#!/usr/bin/env python3
from Bio import SeqIO
import argparse

def convert(input_file, output_file):
    SeqIO.convert(input_file, "genbank", output_file, "fasta")

def get_args():
    parser = argparse.ArgumentParser(description="Convert Genbank files to FASTA format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='Input Genbank file') 
    parser.add_argument('-o', '--output', required=True,help='Output FASTA file')
    return parser.parse_args()

def main():
    args = get_args()
    convert(args.input, args.output)

if __name__ == '__main__':
    main()