#!/usr/bin/env python
"""
Given both reference FASTA files of strain-A and strain-B, a dict mapping
strain-B chromosome names to strain-A chromosome names, and a BED file
containing imprecise structural variants generated by calling
`pbsv run` taking raw reads of strain-A and reference of strain-B,
zoon in to substrings that span each imprecise structural variant on both
reference FASTA files and call `pbsv run` to make a precise structural variant.

e.g., given yeast fy1679 reference FASTA, yeast R64 reference FASTA,
a dict of R64 reference names mapping to fy1679 reference names,
and a bed file of imprecise structural variants from `pbsv run`,
zoom in to substrings spanning each structural variant on both fy1679,
and R64 sequences, then use `pbsv run` to call structural variants.
"""

from argparse import ArgumentParser
import os.path as op
import sys
import os
import logging
from collections import defaultdict
from pbcore.io import DataSet, FastaWriter

from pbsv.independent.utils import execute, realpath, execute_as_bash
from pbsv.ngmlrmap import make_fai_cmd
from pbsv.cli import _mkdir
from pbsv.io.VcfIO import BedReader

from .independent import Constants as C
from .io import SVPolishFiles
from .utils import write_to_bash_file,get_query_subreads_from_alns, get_query_zmws_from_alns, bed2prefix, substr_fasta, basename_prefix_of_fn, get_ref_extension_for_sv, make_subreads_bam_of_zmws2, pbsv_run_and_transform_cmds
from pbsv.libs import Fastafile
from .sv_pbdagcon import get_region_of_seq_in_a_match_b

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()


def polish_a_sv(bed_record, alns, out_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa, make_subreads_bam, make_scripts, execute_scripts, min_qv, ref_ext_len, use_sge):
    """
    Given a structural variant (as bed_record) and its supportive alignments, polish the given structural variant.
    * if make_reference_fa is True, generate a substr of chromosome as reference for this structural variant
    * if make_subreads_bam is True, generate a subreads bam of all subreads of zmws spanning this strucrtural variant
    * if make_scripts is True, generate scripts to call polished structural variant
    """
    srs = get_query_subreads_from_alns(alns)
    zmws = get_query_zmws_from_alns(alns)

    sv_prefix = bed2prefix(bed_record)
    data_dir = realpath(op.join(out_dir, sv_prefix))
    svp_files_obj = SVPolishFiles(root_dir=data_dir, min_qv=min_qv, ref_ext_len=ref_ext_len)
    _mkdir(svp_files_obj.root_dir) # make a subdirectory (e.g., chrI_0_100_Deletion_-100) for all polishing files

    if make_reference_fa:
        # make a substring spanning the expected structural variants
        #ref_start, ref_end = max(0, bed_record.start - C.REFERENCE_EXTENSION), bed_record.end + C.REFERENCE_EXTENSION
        ref_start, ref_end = get_ref_extension_for_sv(bed_record)
        substr_fasta(fileobj=reference_fasta_obj, chrom=bed_record.chrom, start=ref_start, end=ref_end, out_fa_fn=svp_files_obj.sv_ref_fa)

    if make_subreads_bam:
        # get all raw subreads spanning the structural variants
        make_subreads_bam_of_zmws2(in_subreads_fn_or_obj=subreads_ds_obj, zmws=zmws, out_bam_fn=svp_files_obj.subreads_bam, out_fa_fn=svp_files_obj.subreads_fa)

    if make_scripts:
        svp_files_obj.make_all_scripts()

    if execute_scripts:
        svp_files_obj.execute_all_scripts(use_sge=use_sge)

from StringIO import StringIO

def get_dict_from_json(fp):
    """Get a dict from json"""
    import json
    try:
        return dict(json.load(fp))
    except Exception as e:
        raise ValueError("Could not get a dict from %s" % fp.filename)

def make_precise_sv(ref_a_obj, ref_b_obj, bed_record, b2a_names, work_dir):
    """
    ref_a_obj --- Fastafile object of reference FASTA file a
    ref_b_obj --- Fastafile object of reference FASTA file b
    bed_record --- structural variant called by mapping raw reads of a to reference of b
    b2a_names --- a dict mapping sequence names in b to sequence names in a
    """
    _mkdir(work_dir)
    a_chrom, b_chrom = b2a_names[bed_record.chrom], bed_record.chrom
    ref_b_start, ref_b_end = get_ref_extension_for_sv(bed_record, ref_seq_len=ref_b_obj.lengths[ref_b_obj.references.index(b_chrom)])

    # make a FASTA file substr_b_fn containing exactly one sequence of reference b which spans structural variant
    substr_b_fn = op.join(work_dir, '%s.%s_%s_%s.fasta' % (basename_prefix_of_fn(ref_b_obj.filename), b_chrom, ref_b_start, ref_b_end))
    substr_fasta(fileobj=ref_b_obj, chrom=b_chrom, start=ref_b_start, end=ref_b_end, out_fa_fn=substr_b_fn)

    # make a FASTA file substr_a_fn containing exactly one sequence of reference a which matches the sequence in substr_b_fn
    chrom, ref_a_start, ref_a_end = get_region_of_seq_in_a_match_b(a_fa_obj=ref_a_obj, b_fa_obj=Fastafile(substr_b_fn), a_seq_name=a_chrom, work_dir=work_dir)
    assert chrom == a_chrom
    substr_a_fn = op.join(work_dir, '%s.%s_%s_%s.fasta' % (basename_prefix_of_fn(ref_a_obj.filename), a_chrom, ref_a_start, ref_a_end))
    if ref_a_start is not None and ref_a_end is not None:
        substr_fasta(fileobj=ref_a_obj, chrom=a_chrom, start=ref_a_start, end=ref_a_end, out_fa_fn=substr_a_fn)
    else:
        raise ValueError("Could not find any substring in %s that matches the sequence in %s" % (ref_a_obj.filename, ref_b_obj.filename))

    # call `pbsv run` to call structural variants
    o_prefix = op.join(work_dir, '%s.%s' % (basename_prefix_of_fn(substr_a_fn), basename_prefix_of_fn(substr_b_fn)))
    o_bam_fn, o_bed_fn, o_sh_fn = o_prefix + '.chained.bam', o_prefix + '.bed', o_prefix + '.pbsv_run.sh'
    cmds = pbsv_run_and_transform_cmds(reads_fn=substr_a_fn, ref_fa_fn=substr_b_fn, cfg_fn=C.PBSV_POLISH_CFG, o_bam_fn=o_bam_fn, o_bed_fn=o_bed_fn, algorithm='ngmlr')
    write_to_bash_file(cmds=cmds, bash_sh_fn=o_sh_fn)


def run(args):
    ref_a_fn, ref_b_fn, in_bed_fn, b2a_names_json_fn, out_dir, = args.ref_a_fn, args.ref_b_fn, args.in_bed_fn, args.b2a_names_json_fn, args.out_dir
    _mkdir(out_dir)

    # bed_objs = get_sv_from_non_pbsv_bed(in_bed_fn)
    bed_reader = BedReader(in_bed_fn)
    ref_a_obj = Fastafile(ref_a_fn)
    ref_b_obj = Fastafile(ref_b_fn)
    b2a_names = get_dict_from_json(open(b2a_names_json_fn, 'r'))

    for bed_record in bed_reader:
        sv_prefix = bed2prefix(bed_record)
        data_dir = realpath(op.join(out_dir, sv_prefix))
        log.info("Processing sv %s" % sv_prefix)
        print("Processing sv %s" % sv_prefix)
        make_precise_sv(ref_a_obj, ref_b_obj, bed_record, b2a_names, data_dir)

    bed_reader.close()
    ref_a_obj.close()
    ref_b_obj.close()


def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("")
    parser.add_argument("ref_a_fn", type=str, help="Reference FASTA file of strain-A")
    parser.add_argument("ref_b_fn", type=str, help="Reference FASTA file of strain-B")
    parser.add_argument("in_bed_fn", type=str, help="BED file containing imprecise structural variants mapping strain-A raw reads to strain-B reference")
    parser.add_argument("b2a_names_json_fn", type=str, help="A JSON file containing a dict mapping strain-B chromosome names to strain-B chromosome names")
    parser.add_argument("out_dir", type=str, help="Output Directory")
    #parser.add_argument("--use_sge", default=False, action='store_true', help="If True, use SGE; otherwise, run locally")
    return parser


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()