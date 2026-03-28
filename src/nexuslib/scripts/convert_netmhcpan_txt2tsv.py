"""
The purpose of this python3 file is to convert NetMHCpan TXT output file to a TSV file.
"""


import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Convert NetMHCpan TXT output file to a TSV file.")
    parser.add_argument("-i", dest='input_txt_file', help="NetMHCpan TXT output file.")
    parser.add_argument("-o", dest='output_tsv_file',  help="Output TSV file.")
    arguments = parser.parse_args()
    return arguments


def convert_netmhcpan_txt_to_dataframe(txt_file) -> pd.DataFrame:
    # Step 1. Capture headers
    headers = []
    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith("Pos") and "core" in line:
                headers = line.split()

    # Step 2. Capture data
    data = {}
    for header in headers:
        data[header] = []
    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith("Pos") and "core" in line:
                continue
            if line.startswith("---") or line.startswith("Protein"):
                continue
            if 'Distance to training data' in line:
                continue
            elements = line.split()
            if len(elements) == 16:
                for i,element in enumerate(elements):
                    data[headers[i]].append(element)
                data[headers[-1]].append('')
            elif len(elements) == 18:
                for i in range(0,16):
                    data[headers[i]].append(elements[i])
                data[headers[-1]].append(elements[-1])
            else:
                raise Exception('Unsupported number of data elements: %i: %s' % (len(elements), elements))
    df = pd.DataFrame(data)
    return df


def run():
    args = parse_args()
    df = convert_netmhcpan_txt_to_dataframe(txt_file=args.input_txt_file)
    df["Aff(nM)"] = df["Aff(nM)"].astype(float)
    df.to_csv(args.output_tsv_file, sep='\t', index=False)