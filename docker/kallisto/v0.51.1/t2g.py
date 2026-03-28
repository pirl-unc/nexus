#!/usr/bin/env python

import sys, argparse

def create_transcript_list(input_file, use_name=True, use_version=True):
    r = {}
    for line in input_file:
        if len(line) == 0 or line[0] == '#':
            continue
        l = line.strip().split('\t')
        if l[2] == 'transcript':
            info = l[8]
            d = {}
            for x in info.split('; '):
                x = x.strip()
                p = x.find(' ')
                if p == -1:
                    continue
                k = x[:p]
                p = x.find('"', p)
                p2 = x.find('"', p + 1)
                v = x[p + 1:p2]
                d[k] = v

            if 'transcript_id' not in d or 'gene_id' not in d:
                continue

            tid = d['transcript_id'].split(".")[0]
            gid = d['gene_id'].split(".")[0]
            if use_version:
                if 'transcript_version' not in d or 'gene_version' not in d:
                    continue

                # Uncomment if you want to add versions
                # tid += '.' + d['transcript_version']
                # gid += '.' + d['gene_version']
            gname = None
            if use_name:
                if 'gene_name' not in d:
                    continue
                gname = d['gene_name']

            if tid in r:
                continue

            r[tid] = (gid, gname)
    return r


def print_output(output, r, use_name=True):
    for tid in r:
        if use_name:
            output.write(f"{tid}\t{r[tid][0]}\t{r[tid][1]}\n")
        else:
            output.write(f"{tid}\t{r[tid][0]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates transcript-to-gene info from GTF files.')
    parser.add_argument('gtf_file', type=argparse.FileType('r'), help='Input GTF file')
    parser.add_argument('--use_version', '-v', action='store_true', help='Use version numbers in transcript and gene IDs')
    parser.add_argument('--skip_gene_names', '-s', action='store_true', help='Do not output gene names')
    args = parser.parse_args()

    r = create_transcript_list(args.gtf_file, use_name=not args.skip_gene_names, use_version=args.use_version)
    print_output(sys.stdout, r)
