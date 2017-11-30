#!/usr/bin/env python

from argparse import ArgumentParser
import sys

from pbsv.independent.utils import is_bed, is_vcf
from pbsv.io.VcfIO import VcfReader, VcfWriter, VcfRecord, BedReader, BedWriter, BedRecord


def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser()
    parser.add_argument("input_sv_fn", help="Input BED or VCF filename")
    parser.add_argument("output_sv_fn", help="Output BED or VCF filename, format must be the same as input")
    return parser


def get_reader(fn):
    """Return a BedReader obj or VcfReader obj depending on input fn format"""
    if is_bed(fn):
        return BedReader(fn)
    elif is_vcf(fn):
        return VcfReader(fn)
    else:
        raise ValueError("Could not get reader for %s" % fn)


def get_writer(fn, samples):
    """Return a BedWriter obj or VcfWriter obj depending on fn format"""
    if is_bed(fn):
        return BedWriter(fn, samples)
    elif is_vcf(fn):
        return VcfWriter(fn, samples)
    else:
        raise ValueError("Could not get reader for %s" % fn)


def transform_coordinate_of_sv(svobj):
    """
    Given an svobj, which must either be BedRecord or a VcfRecord, where its chrom must
    conform to format `chrom__substr__start_end`, (e.g., `chr01__substr__11838__13838`),
    transform its coordinate to `chrom`.
    ...doctest:
        >>> from pbsv.independent.common import SvFmts
        >>> o1 = BedRecord(chrom='chr1__substr__100_200', start=0, end=100, sv_type='Deletion', sv_id=None, \
                sv_len=-100, alt=None, fmts=SvFmts.fromDict({'SAMPLE1': '0/1:3:6'}), annotations=['ALU'])
        >>> o2 = transform_coordinate_of_sv(o1)
        >>> o2.chrom, o2.start, o2.end
        ('chr1', 100, 200)
    """
    chrom, start, end = get_chrom_start_end_from_string(svobj.chrom)
    if isinstance(svobj, BedRecord):
        return BedRecord(chrom=chrom, start=svobj.start + start, end=svobj.end + start,
                         sv_type=svobj.sv_type, sv_id=svobj.sv_id,
                         sv_len=svobj.sv_len, alt=svobj.alt, annotations=svobj.annotations,
                         fmts=svobj.fmts)
    elif isinstance(svobj, VcfRecord):
        return VcfRecord(chrom=chrom, pos=svobj.pos + start, end=svobj.end + start,
                         ref=svobj.ref, alt=svobj.alt, fmts=svobj.fmts,
                         annotations=svobj.annotations, sv_type=svobj.sv_type, sv_len=svobj.sv_len)
    else:
        raise TypeError("svobj must be either BedRecord or VcfRecord while it is {}".format(type(svobj)))


def get_chrom_start_end_from_string(s):
    """Get chrom name, int(start), int(end) from a string '{chrom}__substr__{start}_{end}'
    ...doctest:
    >>> get_chrom_start_end_from_string('chr01__substr__11838_13838')
    ('chr01', 11838, 13838)
    """
    try:
        chrom, s_e = s.split('__substr__')
        start, end = s_e.split('_')
        return chrom, int(start), int(end)
    except Exception:
        raise ValueError("String %s must be of format '{chrom}__substr__{start}_{end}'" % s)


def transform_coordinates(i_fn, o_fn):
    """Transform coordinates of all structural variants in i_fn to o_fn"""
    if not (all(is_bed(fn) for fn in [i_fn, o_fn]) or all(is_vcf(fn) for fn in [i_fn, o_fn])):
        raise ValueError("Input and output must be both BED or VCF")
    with get_reader(i_fn) as reader:
        with get_writer(o_fn, reader.samples) as writer:
            for r in reader:
                writer.writeRecord(transform_coordinate_of_sv(r))


def run(args):
    i_fn, o_fn = args.input_sv_fn, args.output_sv_fn
    transform_coordinates(i_fn, o_fn)


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))


if __name__ == "__main__":
    main()
