# """
# The purpose of this python3 file is to add 'chr' prefix to  an ABRA2 targets BED file.
#
# Last updated date: July 8, 2022
#
# Author: Jin Seok (Andy) Lee
# """
#
#
# import sys
#
#
# def add_chr_to_gtf(input_file, output_file):
#     """
#     Add 'chr' prefix to specific chromosome names (1-22, X, Y) in a GTF file
#     """
#     # Define the chromosomes we want to modify
#     target_chromosomes = set([str(i) for i in range(1, 23)] + ['X', 'Y'])
#
#     with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
#         for line in infile:
#             # Skip header lines
#             if line.startswith('#'):
#                 outfile.write(line)
#                 continue
#
#             # Split the line into fields
#             fields = line.strip().split('\t')
#
#             # Check if we have enough fields (GTF should have 9 fields)
#             if len(fields) >= 9:
#                 chromosome = fields[0]
#
#                 # Add 'chr' only to target chromosomes that don't already have it
#                 if chromosome in target_chromosomes and not chromosome.startswith('chr'):
#                     fields[0] = 'chr' + chromosome
#
#                 # Write the modified line
#                 outfile.write('\t'.join(fields) + '\n')
#             else:
#                 # Write malformed lines as-is
#                 outfile.write(line)
#
#
# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python add_chr_to_gtf_file.py input.gtf output.gtf")
#         sys.exit(1)
#
#     input_file = sys.argv[1]
#     output_file = sys.argv[2]
#
#     add_chr_to_gtf(input_file, output_file)
#     print(f"Modified GTF written to {output_file}")
#     print("Added 'chr' prefix to chromosomes: 1-22, X, Y only")
