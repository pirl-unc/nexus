"""
Run isONform separately on each cluster of reads and consolidate the results.

Takes a FASTQ file plus a TSV assigning reads to clusters, assembles each cluster
independently in parallel, and merges the per-cluster isoforms into a single set of
output files.
"""


import argparse
import ast
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
# the parent (Step 5) and read by pool workers WITHOUT being pickled per task.
READS: Dict[str, Tuple[str, str]] = {}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--isonform-parallel",
        dest="isonform_parallel",
        default="isONform_parallel",
        type=str,
        help="Path to (or name on PATH of) the isONform_parallel executable "
             "(default: 'isONform_parallel')."
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
        default="isonform",
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
        "--num-parallel",
        dest="num_parallel",
        type=int,
        default=64,
        help="Number of clusters to assemble concurrently (default: 64). isONform "
             "assembles a single cluster single-threaded, so this is the ONLY "
             "parallelism knob; size the node for roughly num_parallel CPUs."
    )
    parser.add_argument(
        "--extra-args",
        dest="extra_args",
        default='',
        type=str,
        help="Extra isONform_parallel arguments, e.g. \"--k 20 --w 31 --xmin 14 "
             "--xmax 80 --exact_instance_limit 50 --max_seqs_to_spoa 200 --delta_len 10 "
             "--iso_abundance 1 --delta_iso_len_3 30 --delta_iso_len_5 50\" (default: ''). "
             "Do NOT pass --write_fastq: this script reads back transcriptome.fasta."
    )
    args = parser.parse_args()
    return args


def check_dependencies(isonform_parallel: str):
    missing = [
        name for name in (isonform_parallel, "spoa", "python")
        if shutil.which(name) is None
    ]
    if missing:
        raise FileNotFoundError(
            "Not found on PATH: %s. isONform requires the isONform_parallel "
            "executable, the `spoa` binary, and a `python` executable."
            % ", ".join(missing))


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


def make_transcript_id(
        cluster_id: str,
        isonform_id: str
) -> str:
    """
    Rewrite an isONform sequence name into a globally unique transcript ID.

    isONform names its isoforms '<cluster>_<batch>_<consensus>', where <cluster> is
    parsed from the input FASTQ's basename. Every cluster is handed its own input
    directory containing a single '0.fastq' (see run_isonform), so <cluster> is always
    '0' and carries no information; we drop it and substitute the real cluster ID.
    <batch> is retained because large clusters are split into several batches and it
    is what keeps their isoform names distinct.
    """
    parts = isonform_id.split("_", 1)
    suffix = parts[1] if len(parts) == 2 else isonform_id
    return "cid_%s_%s" % (cluster_id, suffix)


def format_isonform_mapping_file(
        mapping_file: Path,
        cluster_id: str
) -> str:
    """
    Read isONform's transcriptome_mapping.txt and return the read-support rows as a
    TSV text blob (no header), one line per (transcript, read) pair, columns in
    READ_SUPPORT_COLUMNS order.

    The file is NOT tabular. It is a two-line-per-record, FASTA-like format whose
    second line is a Python list repr of the supporting read names:

        >0_0_3
        ['read_a', 'read_b', 'read_c']

    hence ast.literal_eval (which parses literals only and never executes code).

    Returning preformatted text (rather than a pandas DataFrame) keeps pandas out of
    the up-to-100k-iteration consolidation loop on the parent, and shrinks what is
    pickled back from each worker.
    """
    lines: List[str] = []
    transcript_id = None
    with open(mapping_file, "rt") as file:
        for line in file:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                transcript_id = make_transcript_id(cluster_id=cluster_id, isonform_id=line[1:])
                continue
            if transcript_id is None:
                continue  # payload line with no preceding header; nothing to attach it to
            read_names = ast.literal_eval(line)
            num_supporting_reads = str(len(read_names))
            for read_name in read_names:
                lines.append("\t".join([
                    cluster_id,             # cluster_id
                    str(read_name),         # read_name
                    transcript_id,          # transcript_id
                    num_supporting_reads,   # num_supporting_reads
                ]))
            transcript_id = None
    if not lines:
        return ""
    return "\n".join(lines) + "\n"


def write_fastq_file(
        reads: Dict[str, Tuple[str, str]],
        fastq_file: Path
):
    with open(fastq_file, "wt") as fh:
        for read_name, (sequence, quality) in reads.items():
            fh.write(f"@{read_name}\n{sequence}\n+\n{quality}\n")


def run_isonform(
        cluster_id: str,
        reads: Dict[str, Tuple[str, str]],
        temp_dir: Path,
        isonform_parallel: str,
        extra_args: str
) -> Tuple[str, bool, str, str]:
    """
    Assemble a single cluster with isONform in an isolated temp directory.

    Returns:
        Tuple[cluster_id, success, fasta_text, mapping_text]
            fasta_text:   ">cid_<cluster>_<batch>_<n>\\nseq\\n" blob to append ("" on failure).
            mapping_text: read-support rows as a TSV text blob ("" on failure).
    """
    with tempfile.TemporaryDirectory(dir=temp_dir, prefix=f"isonform_{cluster_id}") as tmp:
        tmpdir = Path(tmp)
        temp_input_dir = tmpdir / "in"
        temp_input_dir.mkdir(parents=True, exist_ok=True)
        write_fastq_file(reads=reads, fastq_file=temp_input_dir / "0.fastq")

        temp_output_dir = tmpdir / "out"
        temp_output_dir.mkdir(parents=True, exist_ok=True)
        temp_split_dir = tmpdir / "split"
        temp_split_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            isonform_parallel,
            "--fastq_folder", str(temp_input_dir),
            "--outfolder", str(temp_output_dir),
            "--tmpdir", str(temp_split_dir),
            "--t", "1",
            "--split_wrt_batches",
            *shlex.split(extra_args)] # splits tokens; empty -> nothing

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print('cluster ID %s: isONform exit %s' % (cluster_id, proc.returncode))
            print(proc.stdout)
            print(proc.stderr)
            return (cluster_id, False, "", "")

        transcripts_file = temp_output_dir / "transcriptome.fasta"
        if not transcripts_file.exists():
            print("cluster ID %s produced no transcriptome.fasta file." % cluster_id)
            return (cluster_id, False, "", "")

        mapping_file = temp_output_dir / "transcriptome_mapping.txt"
        if not mapping_file.exists():
            print("cluster ID %s produced no transcriptome_mapping.txt file." % cluster_id)
            return (cluster_id, False, "", "")

        fasta_lines = []
        with pysam.FastxFile(str(transcripts_file)) as fh_in:
            for entry in fh_in:
                transcript_id = make_transcript_id(
                    cluster_id=cluster_id,
                    isonform_id=str(entry.name)
                )
                fasta_lines.append(f">{transcript_id}\n{entry.sequence}\n")
        fasta_text = "".join(fasta_lines)

        mapping_text = format_isonform_mapping_file(
            mapping_file=mapping_file,
            cluster_id=cluster_id
        )

    return (cluster_id, True, fasta_text, mapping_text)


def run_isonform_task(
        task: Tuple[str, List[str], Path, str, str]
) -> Tuple[str, bool, str, str]:
    cluster_id, read_names, temp_dir, isonform_parallel, extra_args = task
    try:
        reads = {read_name: READS[read_name] for read_name in read_names}
        return run_isonform(
            cluster_id=cluster_id,
            reads=reads,
            temp_dir=temp_dir,
            isonform_parallel=isonform_parallel,
            extra_args=extra_args
        )
    except Exception as e:
        print('cluster ID %s failed: %s' % (cluster_id, e))
        return (cluster_id, False, "", "")


if __name__ == "__main__":
    # Step 1. Parse input arguments.
    args = parse_args()

    # Step 2. Verify isONform and its external dependencies are runnable.
    check_dependencies(isonform_parallel=args.isonform_parallel)

    # Step 3. Create output directory.
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Step 4. Resolve temp directory.
    temp_dir = resolve_temp_dir(temp_dir=args.temp_dir)
    print("Temp dir: %s" % temp_dir)

    # Step 5. Resolve output file names.
    if args.output_prefix == "":
        args.output_prefix = "isonform"
    else:
        args.output_prefix = args.output_prefix + "_isonform"

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
                args.isonform_parallel,
                args.extra_args
            )

    n = len(ordered_cluster_ids)
    num_parallel = max(1, args.num_parallel)
    cpu_count = os.cpu_count() or 1
    if num_parallel > cpu_count:
        print("WARNING: num_parallel(%i) exceeds host CPUs (%i); expect oversubscription." %
              (num_parallel, cpu_count))
    elif num_parallel < cpu_count:
        print("NOTE: num_parallel(%i) leaves host CPUs (%i) idle; isONform assembles "
              "each cluster single-threaded, so raising --num-parallel is the only way "
              "to use them." % (num_parallel, cpu_count))
    print("Assembling %i clusters with %i parallel worker(s)." % (n, num_parallel))

    # Step 8. Assemble each cluster, consolidating results in the parent.
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
            results = map(run_isonform_task, iter_tasks()) # serial; no Pool
        else:
            ctx = multiprocessing.get_context("fork")
            pool = ctx.Pool(processes=num_parallel, maxtasksperchild=100)
            results = pool.imap_unordered(run_isonform_task, iter_tasks())

        for cluster_id, success, fasta_text, mapping_text in results:
            completed += 1
            if completed % 1000 == 0:
                print("[%s] Completed %i/%i clusters." %
                      (datetime.now().strftime("%m/%d/%Y %H:%M:%S"), completed, n))
            if success:
                if fasta_text:
                    fh_fasta.write(fasta_text)
                if mapping_text:
                    fh_tsv.write(mapping_text)
            data['cluster_id'].append(cluster_id)
            data['status'].append("succeeded" if success else "failed")

    if pool is not None:
        pool.close()
        pool.join()

    # Step 9. Write the per-cluster status table.
    output_status_tsv_file = args.output_dir / ("%s_status.tsv" % args.output_prefix)
    df_status = pd.DataFrame(data)
    df_status.to_csv(output_status_tsv_file, sep="\t", index=False)