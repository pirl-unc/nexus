import argparse
import pandas as pd
import yaml
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description="Create Beers2 input data."
    )
    ap.add_argument(
        "--sample-id",
        required=True,
        type=str,
        dest="sample_id",
        help="Sample ID."
    )
    ap.add_argument(
        "--tsv-file",
        required=True,
        type=str,
        dest="tsv_file",
        help="Input TSV file with the following headers: 'transcript_id', 'num_molecules', 'sequence'."
    )
    ap.add_argument(
        "--config-yaml-file",
        required=True,
        type=str,
        dest="config_yaml_file",
        help="BEERS2 config yaml file."
    )
    ap.add_argument(
        "--output-dir",
        required=True,
        type=str,
        dest="output_dir",
        help="Output directory."
    )
    ap.add_argument(
        "--polya-length",
        type=int,
        dest="polya_length",
        default=150,
        help="PolyA tail length to append (0 to disable). BEERS2's default PolyAStep filters out untailed molecules (default: 150)."
    )
    return ap.parse_args()


def run():
    # Step 1. Parse input arguments
    args = parse_args()

    # Step 2. Create necessary output directories
    output_dir = args.output_dir + '/' + args.sample_id + '/'
    Path(output_dir + '/sample1/').mkdir(parents=True, exist_ok=True)

    # Step 3. Create molecules TXT file
    output_txt_file = output_dir + '/sample1/molecule_0.txt'
    df = pd.read_csv(args.tsv_file, sep="\t")
    with open(output_txt_file, "w") as out:
        for _, row in df.iterrows():
            transcript_id = row["transcript_id"]
            num_molecules = row["num_molecules"]
            sequence = row["sequence"]
            sequence = sequence + ("A" * args.polya_length)
            L = len(sequence)
            cigar = f"{L}M"
            for i in range(1, num_molecules + 1):
                # transcript_id needs uniqueness if you emit copies
                molecule_id = f"{transcript_id}_{i}"
                row = [molecule_id, transcript_id, "1", cigar, "1", cigar, "+", sequence]
                out.write("\t".join(row) + "\n")

    # Step 4. Create a dummy reference genome FASTA file
    # @SQ headers for SAM output (chromosome names + lengths). Extract
    # the chromosome names referenced by the molecules and emit a stub
    # sequence of N's long enough to cover the longest mapped molecule.
    chrom_to_max_end = {}
    with open(output_txt_file) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            chrom = cols[1]
            # Molecule sequence is the 8th column; use its length as a
            # safe upper bound for the chromosome end coordinate.
            seq_len = len(cols[7])
            chrom_to_max_end[chrom] = max(chrom_to_max_end.get(chrom, 0), seq_len)
    output_reference_fasta_file = output_dir + 'reference.fa'
    with open(output_reference_fasta_file, 'w') as fh:
        for chrom, length in sorted(chrom_to_max_end.items()):
            # Pad generously to avoid off-by-one issues with @SQ LN.
            fh.write('>%s\n%s\n' % (chrom, 'N' * (length + 1000)))

    # Step 5. Create config yaml file
    with open(args.config_yaml_file) as fh:
        cfg = yaml.safe_load(fh)
        cfg['global_config']['resources']['reference_genome_fasta'] = output_reference_fasta_file
        cfg['library_prep_pipeline']['input']['directory_path'] = args.output_dir + '/' + args.sample_id + '/'
    output_config_yaml_file = args.output_dir + '/' + args.sample_id + '/beers2.config.yaml'
    with open(output_config_yaml_file, 'w') as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)

