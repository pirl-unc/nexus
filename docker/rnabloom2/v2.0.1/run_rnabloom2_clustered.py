import argparse
import gzip
import multiprocessing
import os
import pandas as pd
import pysam
import shlex
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


READ_SUPPORT_COLUMNS = [
    'group_id',
    'read_name',
    'read_length',
    'read_start',
    'read_end',
    'strand',
    'transcript_id',
    'transcript_length',
    'transcript_start',
    'transcript_end',
    'num_residue_matches',
    'frac_transcript_covered',
    'alignment_block_length',
    'mapping_quality'
]


# All input reads, keyed by read ID -> (sequence, quality_string). Populated once in
# the parent (Step 5) and read by pool workers WITHOUT being pickled per task.
READS: Dict[str, Tuple[str, str]] = {}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--jar",
        dest="jar",
        required=True,
        type=Path,
        help="Path to RNA-Bloom.jar"
    )
    parser.add_argument(
        "--fastq-file",
        default="fastq_file",
        required=True,
        type=Path,
        help="Input FASTQ file."
    )
    parser.add_argument(
        "--tsv-file",
        default="tsv_file",
        required=True,
        type=Path,
        help="Input TSV file of grouped read IDs. Expected columns: 'cluster_id', 'read_name'."
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
        default="rnabloom2",
        type=str,
        help="Output prefix (default: '')."
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
        help="Threads per RNA-Bloom2 call (i.e. per parallel worker)."
    )
    parser.add_argument(
        "--num-parallel",
        dest="num_parallel",
        type=int,
        default=64,
        help="Number of clusters to assemble concurrently (default: 1). Each "
             "worker runs one RNA-Bloom2 process using --num-threads threads and "
             "--xmx heap, so size the node for roughly num_parallel * num_threads "
             "CPUs and num_parallel * xmx memory."
    )
    parser.add_argument(
        "--xmx",
        dest="xmx",
        default="2g",
        type=str,
        help="JVM heap per RNA-Bloom2 call (i.e. per parallel worker)."
    )
    parser.add_argument(
        "--extra-args",
        dest="extra_args",
        default='',
        type=str,
        help="Extra RNA-Bloom2 arguments (default: '')."
    )
    args = parser.parse_args()
    return args


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
        needed      :   If set, only load reads whose ID is in this set.

    Returns:
        Dict[read ID, Tuple[sequence, quality_string]]
    """
    reads: Dict[str, Tuple[str, str]] = {}
    with pysam.FastxFile(str(fastq_file)) as fh:
        for entry in fh:
            name = str(entry.name)
            if needed is not None and name not in needed:
                continue
            reads[name] = (str(entry.sequence), entry.quality)
    return reads


def format_rnabloom2_paf_file(
        paf_file: Path,
        group_id: str
) -> str:
    lines: List[str] = []
    with gzip.open(str(paf_file), 'rt') as file:
        for line in file:
            values = line.rstrip('\n').split('\t')
            frac_transcript_covered = float(values[9]) / float(values[6])
            row = "\t".join([
                group_id,                                       # group_id
                values[0],                                      # read_name
                values[1],                                      # read_length
                values[2],                                      # read_start
                values[3],                                      # read_end
                values[4],                                      # strand
                "gid_%s_%i" % (group_id, int(values[5])),       # transcript_id
                values[6],                                      # transcript_length
                values[7],                                      # transcript_start
                values[8],                                      # transcript_end
                values[9],                                      # num_residue_matches
                str(frac_transcript_covered),                   # frac_transcript_covered
                values[10],                                     # alignment_block_length
                values[11],                                     # mapping_quality
            ])
            lines.append(row)
    if not lines:
        return ""
    return "\n".join(lines) + "\n"


def write_fastq_file(
        reads: Dict[str, Tuple[str, str]],
        fastq_file: Path
):
    with open(fastq_file, "wt") as fh:
        for read_id, (sequence, quality) in reads.items():
            fh.write(f"@{read_id}\n{sequence}\n+\n{quality}\n")


def run_rnabloom2(
        group_id: str,
        reads: Dict[str, Tuple[str, str]],
        temp_dir: Path,
        rnabloom_jar: Path,
        extra_args: str,
        num_threads: int,
        xmx: str
) -> Tuple[str, bool, str, str]:
    """
    Assemble a single cluster with RNA-Bloom2 in an isolated temp directory.

    Returns:
        Tuple[group_id, success, fasta_text, paf_text]
            fasta_text: ">gid_<group>_<name>\\nseq\\n" blob to append ("" on failure).
            paf_text:   read-support rows as a TSV text blob ("" on failure).
    """
    with tempfile.TemporaryDirectory(dir=temp_dir, prefix=f"rnabloom2_{group_id}") as tmp:
        tmpdir = Path(tmp)
        temp_fastq_file = tmpdir / "reads.fastq"
        write_fastq_file(reads=reads, fastq_file=temp_fastq_file)
        temp_output_dir = tmpdir / "out"
        temp_output_dir.mkdir(parents=True, exist_ok=True)

        # -XX:ActiveProcessorCount pins this JVM to its own thread budget. Without it,
        # every RNA-Bloom2 process sees ALL host cores and sizes its GC + internal
        # pools off that.
        cmd = [
            "java", f"-Xmx{xmx}", f"-XX:ActiveProcessorCount={num_threads}",
            "-jar", str(rnabloom_jar),
            "-long", str(temp_fastq_file),
            "--threads", str(num_threads),
            "--outdir", str(temp_output_dir),
            *shlex.split(extra_args)]  # splits tokens; empty -> nothing

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print('group ID %s: RNAbloom2 exit %s' % (group_id, proc.returncode))
            print(proc.stdout)
            return (group_id, False, "", "")

        transcripts_file = temp_output_dir / "rnabloom.longreads.assembly4.pol.fa"
        if not transcripts_file.exists():
            print("group ID %s produced no rnabloom.longreads.assembly4.pol.fa file." % group_id)
            return (group_id, False, "", "")

        paf_file = temp_output_dir / "rnabloom.longreads.assembly3.map.paf.gz"
        if not paf_file.exists():
            print("group ID %s produced no rnabloom.longreads.assembly3.map.paf.gz file." % group_id)
            return (group_id, False, "", "")

        # Read transcripts into an in-memory FASTA blob before the temp dir is cleaned.
        fasta_lines = []
        with pysam.FastxFile(str(transcripts_file)) as fh_in:
            for entry in fh_in:
                fasta_lines.append(f">gid_{group_id}_{entry.name}\n{entry.sequence}\n")
        fasta_text = "".join(fasta_lines)

        # Read support for this group, preformatted as a TSV text blob.
        paf_text = format_rnabloom2_paf_file(paf_file=paf_file, group_id=group_id)

    return (group_id, True, fasta_text, paf_text)


def run_rnabloom2_task(
        task: Tuple[str, List[str], Path, Path, str, int, str]
) -> Tuple[str, bool, str, str]:
    group_id, read_ids, temp_dir, rnabloom_jar, extra_args, num_threads, xmx = task
    try:
        reads = {read_id: READS[read_id] for read_id in read_ids}
        return run_rnabloom2(
            group_id=group_id,
            reads=reads,
            temp_dir=temp_dir,
            rnabloom_jar=rnabloom_jar,
            extra_args=extra_args,
            num_threads=num_threads,
            xmx=xmx
        )
    except Exception as e:
        print('group ID %s failed: %s' % (group_id, e))
        return (group_id, False, "", "")


if __name__ == "__main__":
    # Step 1. Parse input arguments.
    args = parse_args()

    # Step 2. Create output directory.
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Step 3. Resolve temp directory.
    temp_dir = resolve_temp_dir(temp_dir=args.temp_dir)
    print("Temp dir: %s" % temp_dir)

    # Step 4. Resolve output file names.
    if args.output_prefix == "":
        args.output_prefix = "rnabloom2"
    else:
        args.output_prefix = args.output_prefix + "_rnabloom2"

    output_fasta_file = args.output_dir / ("%s_merged.fasta.gz" % args.output_prefix)
    output_reads_tsv_file = args.output_dir / ("%s_merged.read_support.tsv" % args.output_prefix)

    # Step 5. Load the groups and the reads.
    # Only reads referenced by some cluster are loaded.
    df_clusters = pd.read_csv(args.tsv_file, sep='\t')
    needed_read_ids = set(df_clusters['read_name'].astype(str).unique())
    READS = load_fastq_file(fastq_file=args.fastq_file, needed=needed_read_ids)

    # Step 6. Order groups largest-first, then build a LAZY task stream in that order.
    grouped_read_ids = df_clusters.groupby('cluster_id')['read_name']
    ordered_group_ids = (
        grouped_read_ids.nunique().sort_values(ascending=False).index.tolist()
    )

    def iter_tasks():
        for group_id in ordered_group_ids:
            read_ids = [str(read_id) for read_id in grouped_read_ids.get_group(group_id).unique()]
            yield (
                str(group_id),
                read_ids,
                temp_dir,
                args.jar,
                args.extra_args,
                args.num_threads,
                args.xmx
            )

    n = len(ordered_group_ids)
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
    print("Assembling %i groups with %i parallel worker(s), %i thread(s) each." %
          (n, num_parallel, args.num_threads))

    # Step 7. Assemble each group, consolidating results in the parent.
    data = {
        'group_id': [],
        'status': []
    }
    completed = 0
    pool = None
    with gzip.open(output_fasta_file, "wt") as fh_fasta, \
            open(output_reads_tsv_file, "w") as fh_tsv:
        fh_tsv.write("\t".join(READ_SUPPORT_COLUMNS) + "\n")

        if num_parallel <= 1:
            results = map(run_rnabloom2_task, iter_tasks())  # serial; no Pool
        else:
            # Force 'fork' so workers inherit READS copy-on-write (see READS).
            ctx = multiprocessing.get_context("fork")
            pool = ctx.Pool(processes=num_parallel, maxtasksperchild=100)
            results = pool.imap_unordered(run_rnabloom2_task, iter_tasks())

        for group_id, success, fasta_text, paf_text in results:
            completed += 1
            if completed % 1000 == 0:
                print("[%s] Completed %i/%i groups." %
                      (datetime.now().strftime("%m/%d/%Y %H:%M:%S"), completed, n))
            if success:
                if fasta_text:
                    fh_fasta.write(fasta_text)
                if paf_text:
                    fh_tsv.write(paf_text)
            data['group_id'].append(group_id)
            data['status'].append("succeeded" if success else "failed")

    if pool is not None:
        pool.close()
        pool.join()

    # Step 8. Write the per-group status table.
    output_status_tsv_file = args.output_dir / ("%s_status.tsv" % args.output_prefix)
    df_status = pd.DataFrame(data)
    df_status.to_csv(output_status_tsv_file, sep="\t", index=False)
