#!/usr/bin/env python

from argparse import ArgumentParser
import sys

from pbsv.independent.utils import is_fasta
from pbsv.libs import Fastafile, AlignmentFile

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("Extract a substring (chrom, start, end) from input FASTA file and save to output FASTA file.")
    parser.add_argument("input_fa_fn", help="Input FASTA filename")
    parser.add_argument("chrom", help="Name of the sequence to extract")
    parser.add_argument("start", help="Start position of substring, 0-based, inclusive")
    parser.add_argument("end", help="End position of substring, 0-based, exclusive")
    parser.add_argument("output_fa_fn", help="Output FASTA filename")
    return parser

from .utils import substr_fasta

def run(args):
    assert is_fasta(args.input_fa_fn) and is_fasta(args.output_fa_fn)
    fileobj = Fastafile(args.input_fa_fn)
    substr_fasta(fileobj=fileobj, chrom=str(args.chrom), start=int(args.start), end=int(args.end), out_fa_fn=args.output_fa_fn)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
