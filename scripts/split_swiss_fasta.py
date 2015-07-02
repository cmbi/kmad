#!/usr/bin/python
import argparse


def split(input_name, output_dir):
    with open(input) as a:
        uni_dat = a.read().splitlines()

    for lineI in uni_dat:
        if ">" == lineI[0]:
            if "out" in locals() and not out.closed:
                out.close()
            uni_id = lineI.split()[0].split("|")[2]
            out = open(output_dir + uni_id + ".fasta", 'w')
        if "out" in locals() and not out.closed:
            out.write(lineI)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a fasta file with"
                                                 "multiple sequences into"
                                                 "single sequence fasta files")
    parser.add_argument('input_name')
    parser.add_argument('output_dir')
    args = parser.parse_args()
    split(args.input_name, args.output_dir)
