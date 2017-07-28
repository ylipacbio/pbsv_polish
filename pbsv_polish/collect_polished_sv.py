#!/usr/bin/env python
from argparse import ArgumentParser
import os
import os.path as op
import sys

from pbsv.independent.utils import execute, realpath
from pbsv.io.VcfIO import BedReader, BedWriter

from .utils import bed2prefix, SVPolishFiles, Constant

def run(args):
    in_bed_fn, out_dir, collected_bed_fn = args.in_bed_fn,  args.out_dir, args.collected_bed_fn
    min_qv = args.min_qv
    writer = BedWriter(collected_bed_fn)

    for bed_record in BedReader(in_bed_fn):
        sv_prefix = bed2prefix(bed_record)
        print("Collecting polished SV for %s" % sv_prefix)
        data_dir = realpath(op.join(out_dir, sv_prefix))
        svp_files_obj = SVPolishFiles(root_dir=data_dir, min_qv=min_qv)
        polished_bed_fn = svp_files_obj.polish_ngmlr_bed

        if not op.exists(polished_bed_fn):
            print("No Polished structural variant detected, use the original one: %s " % ' '.join(str(bed_record).split()[0:5]))
            writer.writeRecord(bed_record)
        else:
            polished_bed_records = [r for r in BedReader(polished_bed_fn)]
            if len(polished_bed_records) == 0:
                print("No Polished structural variant detected, use the original one: %s " % ' '.join(str(bed_record).split()[0:5]))
                writer.writeRecord(bed_record)
            else:
                for r in polished_bed_records:
                    writer.writeRecord(r)

    writer.close()

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("")
    parser.add_argument("in_bed_fn", type=str, help="Structural variants in BED file")
    parser.add_argument("out_dir", type=str, help="Output Directory")
    parser.add_argument("collected_bed_fn", type=str, help="Polished structural variants in BED file")
    parser.add_argument("--min_qv", default=Constant.MIN_POLISH_QV, type=int, help="Minimum Polished QV to include bases in consensus sequence")
    return parser


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
