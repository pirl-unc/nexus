import argparse
import gzip
import multiprocessing
import os
import pandas as pd
import pysam
import shlex
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


READ_SUPPORT_COLUMNS = [
    'cluster_id',
    'read_name',
    'transcript_id',
    'num_supporting_reads'
]


# All input reads, keyed by read name -> (sequence, quality_string). Populated once in
# the parent (Step 6) and read by pool workers WITHOUT being pickled per task.
READS: Dict[str, Tuple[str, str]] = {}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--rattle",
        dest="rattle",
        default="rattle",
        type=str,
        help="Path to (or name on PATH of) the rattle executable (default: 'rattle')."
    )
    parser.add_argument(
        "--fastq-file",
        dest="fastq_file",
        required=True,
        type=Path,
        help="Input FASTQ file."
    )
    parser.add_argument(
        "--tsv-file",
        dest="tsv_file",
        required=True,
        type=Path,
        help="Input TSV file of clustered read names. Expected columns: 'cluster_id', 'read_name'."
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        required=True,
        type=Path,
        help="Output directory."
    )
    parser.add_argument(
        "--output-prefix",
        dest="output_prefix",
        required=True,
        default="rattle",
        type=str,
        help="Output prefix."
    )
    parser.add_argument(
        "--temp-dir",
        dest="temp_dir",
        default=None,
        help="Temp dir (default: $TMPDIR or /tmp)"
    )
    parser.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=1,
        help="Threads per RATTLE call (i.e. per parallel worker)."
    )
    parser.add_argument(
        "--num-parallel",
        dest="num_parallel",
        type=int,
        default=64,
        help="Number of clusters to assemble concurrently (default: 64). Each worker "
             "runs the cluster/correct/polish chain using --num-threads threads, so "
             "size the node for roughly num_parallel * num_threads CPUs."
    )
    parser.add_argument(
        "--cluster-extra-args",
        dest="cluster_extra_args",
        default='',
        type=str,
        help="Extra 'rattle cluster' arguments, e.g. \"--iso --rna\" (default: '')."
    )
    parser.add_argument(
        "--correct-extra-args",
        dest="correct_extra_args",
        default='',
        type=str,
        help="Extra 'rattle correct' arguments, e.g. \"-r 2\" (default: '')."
    )
    parser.add_argument(
        "--polish-extra-args",
        dest="polish_extra_args",
        default='',
        type=str,
        help="Extra 'rattle polish' arguments, e.g. \"--rna\" (default: '')."
    )
    args = parser.parse_args()
    return args


def check_dependencies(rattle: str):
    if shutil.which(rattle) is None:
        raise FileNotFoundError("Not found on PATH: %s" % rattle)


def resolve_temp_dir(temp_dir: Optional[str]) -> Path:
    base = temp_dir or os.environ.get("TMPDIR") or "/tmp"
    p = Path(base)
    if not p.is_dir():
        raise FileNotFoundError(f"Temp directory does not exist: {p}")
    return p


def load_fastq_file(
        fastq_file: Path,
        needed: Optional[set] = None
) -> Dict[str, Tuple[str, str]]:
    """
    Load FASTQ file.

    Parameters:
        fastq_file  :   Path to FASTQ file.
        needed      :   If set, only load reads whose name is in this set.

    Returns:
        Dict[read name, Tuple[sequence, quality_string]]
    """
    reads: Dict[str, Tuple[str, str]] = {}
    with pysam.FastxFile(str(fastq_file)) as fh:
        for entry in fh:
            name = str(entry.name)
            if needed is not None and name not in needed:
                continue
            reads[name] = (str(entry.sequence), entry.quality)
    return reads


def write_fastq_file(
        reads: Dict[str, Tuple[str, str]],
        fastq_file: Path
):
    with open(fastq_file, "wt") as fh:
        for read_name, (sequence, quality) in reads.items():
            fh.write(f"@{read_name}\n{sequence}\n+\n{quality}\n")


def make_transcript_id(
        cluster_id: str,
        cluster_index: str
) -> str:
    """
    Build a globally unique transcript ID. RATTLE numbers its polished transcripts
    from 0 within each invocation, and every cluster gets its own invocation, so the
    input cluster ID must be folded in to keep names unique across the merged output.
    """
    return "cid_%s_%s" % (cluster_id, cluster_index)


def parse_cluster_index(name: str) -> Optional[str]:
    """
    Pull the trailing integer out of one of RATTLE's cluster labels, e.g.
    'transcript_cluster_7' / 'gene_cluster_7' / 'new_cluster_7' / 'cluster_7' -> '7'.
    """
    index = name.rsplit("_", 1)[-1].strip()
    return index if index.isdigit() else None


def parse_polish_summary_file(summary_file: Path) -> Dict[str, str]:
    """
    Read polish_summary.tsv and return {correct-stage cluster index -> polished transcript index}.
    """
    old_to_new: Dict[str, str] = {}
    with open(summary_file, "rt") as file:
        for line in file:
            fields = [field.strip() for field in line.strip().split(",")]
            if len(fields) < 2:
                continue
            old_index = parse_cluster_index(fields[0])
            new_index = parse_cluster_index(fields[-1])
            if old_index is not None and new_index is not None:
                old_to_new[old_index] = new_index
    return old_to_new


def format_rattle_corrected_file(
        corrected_file: Path,
        summary_file: Path,
        cluster_id: str
) -> str:
    old_to_new = parse_polish_summary_file(summary_file=summary_file)

    assignments: List[Tuple[str, str]] = []
    with pysam.FastxFile(str(corrected_file)) as fh:
        for entry in fh:
            fields = str(entry.name).split(",")
            if len(fields) < 2:
                continue  # no cluster annotation appended; nothing to attribute
            cluster_index = parse_cluster_index(fields[-1])
            if cluster_index is None:
                continue
            new_index = old_to_new.get(cluster_index)
            if new_index is None:
                continue  # this consensus did not survive polishing
            assignments.append((fields[0], make_transcript_id(cluster_id, new_index)))

    if not assignments:
        return ""

    num_supporting_reads: Dict[str, int] = {}
    for _, transcript_id in assignments:
        num_supporting_reads[transcript_id] = num_supporting_reads.get(transcript_id, 0) + 1

    lines = [
        "\t".join([
            cluster_id,                                 # cluster_id
            read_name,                                  # read_name
            transcript_id,                              # transcript_id
            str(num_supporting_reads[transcript_id]),   # num_supporting_reads
        ])
        for read_name, transcript_id in assignments
    ]
    return "\n".join(lines) + "\n"


def run_stage(
        cmd: List[str],
        cluster_id: str,
        stage: str
) -> bool:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print('cluster ID %s: rattle %s exit %s' % (cluster_id, stage, proc.returncode))
        print(proc.stdout)
        print(proc.stderr)
        return False
    return True


def run_rattle(
        cluster_id: str,
        reads: Dict[str, Tuple[str, str]],
        temp_dir: Path,
        rattle: str,
        cluster_extra_args: str,
        correct_extra_args: str,
        polish_extra_args: str,
        num_threads: int
) -> Tuple[str, bool, str, str]:
    """
    Run the RATTLE cluster -> correct -> polish chain on a single cluster in an
    isolated temp directory.

    Returns:
        Tuple[cluster_id, success, fasta_text, read_support_text]
            fasta_text:        ">cid_<cluster>_<n>\\nseq\\n" blob to append ("" on failure).
            read_support_text: read-support rows as a TSV text blob ("" on failure).
    """
    with tempfile.TemporaryDirectory(dir=temp_dir, prefix=f"rattle_{cluster_id}") as tmp:
        tmpdir = Path(tmp)
        temp_fastq_file = tmpdir / "reads.fastq"
        write_fastq_file(reads=reads, fastq_file=temp_fastq_file)
        temp_output_dir = tmpdir / "out"
        temp_output_dir.mkdir(parents=True, exist_ok=True)

        # Stage 1. Sub-cluster, to produce the clusters.out that correct needs.
        cmd = [
            rattle, "cluster",
            "--input", str(temp_fastq_file),
            "--output", str(temp_output_dir),
            "-t", str(num_threads),
            *shlex.split(cluster_extra_args)]  # splits tokens; empty -> nothing
        if not run_stage(cmd=cmd, cluster_id=cluster_id, stage="cluster"):
            return (cluster_id, False, "", "")

        clusters_file = temp_output_dir / "clusters.out"
        if not clusters_file.exists():
            print("cluster ID %s produced no clusters.out file." % cluster_id)
            return (cluster_id, False, "", "")

        # Stage 2. Error-correct the reads and build a consensus per sub-cluster.
        cmd = [
            rattle, "correct",
            "--input", str(temp_fastq_file),
            "--clusters", str(clusters_file),
            "--output", str(temp_output_dir),
            "-t", str(num_threads),
            *shlex.split(correct_extra_args)]
        if not run_stage(cmd=cmd, cluster_id=cluster_id, stage="correct"):
            return (cluster_id, False, "", "")

        consensi_file = temp_output_dir / "consensi.fq"
        corrected_file = temp_output_dir / "corrected.fq"
        if not consensi_file.exists() or not corrected_file.exists():
            print("cluster ID %s produced no consensi.fq/corrected.fq file." % cluster_id)
            return (cluster_id, False, "", "")

        if consensi_file.stat().st_size == 0:
            return (cluster_id, True, "", "")

        # Stage 3. Polish the consensi into the final transcripts.
        cmd = [
            rattle, "polish",
            "--input", str(consensi_file),
            "--output-folder", str(temp_output_dir),
            "-t", str(num_threads),
            "--summary",
            *shlex.split(polish_extra_args)]
        if not run_stage(cmd=cmd, cluster_id=cluster_id, stage="polish"):
            return (cluster_id, False, "", "")

        transcriptome_file = temp_output_dir / "transcriptome.fq"
        if not transcriptome_file.exists():
            print("cluster ID %s produced no transcriptome.fq file." % cluster_id)
            return (cluster_id, False, "", "")

        summary_file = temp_output_dir / "polish_summary.tsv"
        if not summary_file.exists():
            print("cluster ID %s produced no polish_summary.tsv file." % cluster_id)
            return (cluster_id, False, "", "")

        # Read transcripts into an in-memory FASTA blob before the temp dir is cleaned.
        fasta_lines = []
        with pysam.FastxFile(str(transcriptome_file)) as fh_in:
            for entry in fh_in:
                cluster_index = parse_cluster_index(str(entry.name))
                if cluster_index is None:
                    continue
                transcript_id = make_transcript_id(cluster_id, cluster_index)
                fasta_lines.append(f">{transcript_id}\n{entry.sequence}\n")
        fasta_text = "".join(fasta_lines)

        # Read support for this cluster, preformatted as a TSV text blob.
        read_support_text = format_rattle_corrected_file(
            corrected_file=corrected_file,
            summary_file=summary_file,
            cluster_id=cluster_id
        )

    return (cluster_id, True, fasta_text, read_support_text)


def run_rattle_task(
        task: Tuple[str, List[str], Path, str, str, str, str, int]
) -> Tuple[str, bool, str, str]:
    (cluster_id, read_names, temp_dir, rattle, cluster_extra_args,
     correct_extra_args, polish_extra_args, num_threads) = task
    try:
        reads = {read_name: READS[read_name] for read_name in read_names}
        return run_rattle(
            cluster_id=cluster_id,
            reads=reads,
            temp_dir=temp_dir,
            rattle=rattle,
            cluster_extra_args=cluster_extra_args,
            correct_extra_args=correct_extra_args,
            polish_extra_args=polish_extra_args,
            num_threads=num_threads
        )
    except Exception as e:
        print('cluster ID %s failed: %s' % (cluster_id, e))
        return (cluster_id, False, "", "")


if __name__ == "__main__":
    # Step 1. Parse input arguments.
    args = parse_args()

    # Step 2. Verify RATTLE is runnable.
    check_dependencies(rattle=args.rattle)

    # Step 3. Create output directory.
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Step 4. Resolve temp directory.
    temp_dir = resolve_temp_dir(temp_dir=args.temp_dir)
    print("Temp dir: %s" % temp_dir)

    # Step 5. Resolve output file names.
    if args.output_prefix == "":
        args.output_prefix = "rattle"
    else:
        args.output_prefix = args.output_prefix + "_rattle"

    output_fasta_file = args.output_dir / ("%s_merged.fasta.gz" % args.output_prefix)
    output_reads_tsv_file = args.output_dir / ("%s_merged.read_support.tsv" % args.output_prefix)

    # Step 6. Load the clusters and the reads.
    # Only reads referenced by some cluster are loaded.
    df_clusters = pd.read_csv(args.tsv_file, sep='\t')
    needed_read_names = set(df_clusters['read_name'].astype(str).unique())
    READS = load_fastq_file(fastq_file=args.fastq_file, needed=needed_read_names)

    # Step 7. Order clusters largest-first, then build a LAZY task stream in that order.
    grouped_read_names = df_clusters.groupby('cluster_id')['read_name']
    ordered_cluster_ids = (
        grouped_read_names.nunique().sort_values(ascending=False).index.tolist()
    )

    def iter_tasks():
        for cluster_id in ordered_cluster_ids:
            read_names = [str(read_name) for read_name in grouped_read_names.get_group(cluster_id).unique()]
            yield (
                str(cluster_id),
                read_names,
                temp_dir,
                args.rattle,
                args.cluster_extra_args,
                args.correct_extra_args,
                args.polish_extra_args,
                args.num_threads
            )

    n = len(ordered_cluster_ids)
    num_parallel = max(1, args.num_parallel)
    cpu_count = os.cpu_count() or 1
    if num_parallel * args.num_threads > cpu_count:
        print("WARNING: num_parallel(%i) * num_threads(%i) = %i exceeds host CPUs (%i); "
              "expect oversubscription." %
              (num_parallel, args.num_threads, num_parallel * args.num_threads, cpu_count))
    elif num_parallel * args.num_threads < cpu_count:
        print("NOTE: num_parallel(%i) * num_threads(%i) = %i leaves host CPUs (%i) idle; "
              "for many small clusters, raising --num-parallel (and lowering "
              "--num-threads) usually improves throughput." %
              (num_parallel, args.num_threads, num_parallel * args.num_threads, cpu_count))
    print("Assembling %i clusters with %i parallel worker(s), %i thread(s) each." %
          (n, num_parallel, args.num_threads))

    # Step 8. Process each cluster, consolidating results in the parent.
    data = {
        'cluster_id': [],
        'status': []
    }
    completed = 0
    pool = None
    with gzip.open(output_fasta_file, "wt") as fh_fasta, \
            open(output_reads_tsv_file, "w") as fh_tsv:
        fh_tsv.write("\t".join(READ_SUPPORT_COLUMNS) + "\n")

        if num_parallel <= 1:
            results = map(run_rattle_task, iter_tasks())  # serial; no Pool
        else:
            # Force 'fork' so workers inherit READS copy-on-write.
            ctx = multiprocessing.get_context("fork")
            pool = ctx.Pool(processes=num_parallel, maxtasksperchild=100)
            results = pool.imap_unordered(run_rattle_task, iter_tasks())

        for cluster_id, success, fasta_text, read_support_text in results:
            completed += 1
            if completed % 1000 == 0:
                print("[%s] Completed %i/%i clusters." %
                      (datetime.now().strftime("%m/%d/%Y %H:%M:%S"), completed, n))
            if success:
                if fasta_text:
                    fh_fasta.write(fasta_text)
                if read_support_text:
                    fh_tsv.write(read_support_text)
            data['cluster_id'].append(cluster_id)
            data['status'].append("succeeded" if success else "failed")

    if pool is not None:
        pool.close()
        pool.join()

    # Step 9. Write the per-cluster status table.
    output_status_tsv_file = args.output_dir / ("%s_status.tsv" % args.output_prefix)
    df_status = pd.DataFrame(data)
    df_status.to_csv(output_status_tsv_file, sep="\t", index=False)
