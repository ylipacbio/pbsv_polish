#!/usr/bin/env python

from argparse import ArgumentParser
import os.path as op
import sys
import os
import logging
from collections import defaultdict
from pbcore.io import DataSet, FastaWriter

from pbsv.independent.utils import execute, realpath, write_to_bash_file, execute_as_bash
from pbsv.ngmlrmap import make_fai_cmd
from pbsv.cli import _mkdir

from .utils import *

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()


def make_diagnose_script_for_pbsv_run(o_dir):
    diagnose_fn = op.join(o_dir, 'diagnose.sh')
    cmd = """
fastalen polished.fasta
fastalen polished.hqlq.fasta
blasr polished.fasta sv_reference_w_extension.fasta --header --maxMatch 15 -m 4
blasr polished.hqlq.fasta sv_reference_w_extension.fasta --header --maxMatch 15 -m 4
"""
    with open(diagnose_fn, 'w') as w:
        w.write(cmd)

def polish_a_sv(bed_record, alns, work_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa, make_subreads_bam, make_scripts, execute_scripts, min_qv, use_sge):
    """
    Given a structural variant (as bed_record) and its supportive alignments, polish the given structural variant.
    * if make_reference_fa is True, generate a substr of chromosome as reference for this structural variant
    * if make_subreads_bam is True, generate a subreads bam of all subreads of zmws spanning this strucrtural variant
    * if make_scripts is True, generate scripts to call polished structural variant
    """
    srs = get_query_subreads_from_alns(alns)
    zmws = get_query_zmws_from_alns(alns)
    svp_files_obj = SVPolishFiles(root_dir=work_dir, min_qv=min_qv)
    _mkdir(svp_files_obj.root_dir) # make a subdirectory (e.g., chrI_0_100_Deletion_-100) for all polishing files

    if make_reference_fa:
        # make a substring spanning the expected structural variants
        #ref_start, ref_end = max(0, bed_record.start - Constant.REFERENCE_EXTENSION), bed_record.end + Constant.REFERENCE_EXTENSION
        ref_start, ref_end = get_ref_extension_for_sv(bed_record, reference_fasta_obj.lengths[reference_fasta_obj.references(bed_record.chrom)])
        substr_fasta(fileobj=reference_fasta_obj, chrom=bed_record.chrom, start=ref_start, end=ref_end, o_fasta_fn=svp_files_obj.sv_ref_fa)

    if make_subreads_bam:
        # get all raw subreads spanning the structural variants
        make_subreads_bam_of_zmws2(in_subreads_fn_or_obj=subreads_ds_obj, zmws=zmws, out_bam_fn=svp_files_obj.subreads_bam, out_fa_fn=svp_files_obj.subreads_fa)

    if make_scripts:
        svp_files_obj.make_all_scripts()

    if execute_scripts:
        svp_files_obj.execute_all_scripts(use_sge=use_sge)


def run(args):
    in_dir, out_dir = args.in_dir, args.out_dir
    _mkdir(out_dir)

    aln_fn = op.join(in_dir, "alignments.bam")
    subreads_xml_fn = op.join(in_dir, "subreads.xml")
    genome_fa = op.join(in_dir, "genome.fa")
    bed_fn = realpath(args.in_bed_fn)
    coverage_fn = op.join(out_dir, 'coverage.txt')

    reference_fasta_obj = Fastafile(genome_fa)
    ofile_obj = open(coverage_fn, 'w')
    alnfile_obj = X2PysamReader(aln_fn)._alignment_file
    subreads_ds_obj = SubreadSet(subreads_xml_fn)
    bedreader_obj = BedReader(bed_fn)

    for bed_record, alns in yield_alns_from_bed_file(alnfile_obj, bedreader_obj=bedreader_obj):
        ofile_obj.write('%s\t%s\n' % (len(get_query_subreads_from_alns(alns)), bed_record)) # write coverage info

        sv_prefix = bed2prefix(bed_record)
        work_dir = realpath(op.join(out_dir, sv_prefix))
        log.info("Processing sv %s" % sv_prefix)
        print("Processing sv %s" % sv_prefix)
        polish_a_sv(bed_record, alns, work_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa=True, make_subreads_bam=True, make_scripts=True, execute_scripts=True, min_qv=args.min_qv, use_sge=args.use_sge)

    reference_fasta_obj.close()
    ofile_obj.close()
    alnfile_obj.close()
    subreads_ds_obj.close()
    bedreader_obj.close()

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("")
    parser.add_argument("in_dir", type=str, help="Input FASTA or FASTQ filename")
    parser.add_argument("in_bed_fn", type=str, help="Structural variants in BED file")
    parser.add_argument("out_dir", type=str, help="Output Directory")
    parser.add_argument("--min_qv", default=Constant.MIN_POLISH_QV, type=int, help="Minimum Polished QV to include bases in consensus sequence")
    parser.add_argument("--use_sge", default=False, action='store_true', help="If True, use SGE; otherwise, run locally")
    return parser


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
