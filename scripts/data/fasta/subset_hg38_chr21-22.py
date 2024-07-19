import os


def extract_chromosomes(fasta_file, output_file, chromosomes):
    with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
        write_mode = False
        for line in infile:
            if line.startswith('>'):
                write_mode = any(chrom in line for chrom in chromosomes)
            if write_mode:
                outfile.write(line)


if __name__ == "__main__":
    fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
    output_fasta_file = '../../../test/data/fasta/hg38_chr21-22.fa'
    extract_chromosomes(fasta_file, output_fasta_file, ['chr21', 'chr22'])
    os.system('gzip ../../../test/data/fasta/hg38_chr21-22.fa')
