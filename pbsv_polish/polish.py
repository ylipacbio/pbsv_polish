#!/usr/bin/env python

import argparse
import os.path as op
import sys
from pbcore.io import FastaWriter

from pbsv.independent.utils import realpath
from pbsv.cli import _mkdir
from pbsv.io.bamstream import SingleFileOpener
from .independent import Constants as C
from .io import SVPolishFiles
from .utils import *

import logging
logging.basicConfig(format='%(asctime) %(message)s')
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()

polish_desc = 'Polish input structural variants'


def polish_a_sv(bed_record, alns, work_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa, make_subreads_bam, make_scripts, execute_scripts, min_qv, ref_ext_len, use_sge, report_f):
    """
    Given a structural variant (as bed_record) and its supportive alignments, polish the given structural variant.
    * if make_reference_fa is True, generate a substr of chromosome as reference for this structural variant
    * if make_subreads_bam is True, generate a subreads bam of all subreads of zmws spanning this strucrtural variant
    * if make_scripts is True, generate scripts to call polished structural variant
    """
    cov = sum([bed_record.fmts[sample].ad for sample in bed_record.samples])
    if cov <= C.MIN_POLISH_COVERAGE:
        msg = "SV={sv}, Coverage={cov}, Skipping structural variant for not having enough coverage!".format(
            cov=cov, sv=bed2prefix(bed_record))
        log.warning(msg)
        report_f(msg)
        return

    # srs = get_query_subreads_from_alns(alns)
    zmws = get_query_zmws_from_alns(alns)
    svp_files_obj = SVPolishFiles(root_dir=work_dir, min_qv=min_qv, ref_ext_len=ref_ext_len)
    _mkdir(svp_files_obj.root_dir)  # make a subdirectory (e.g., chrI_0_100_Deletion_-100) for all polishing files
    # TODO: special treatment for heterzygous sv?
    if make_reference_fa:
        # make a substring spanning the expected structural variants
        ref_start, ref_end = get_ref_extension_for_sv(
            bed_record, reference_fasta_obj.get_reference_length(bed_record.chrom))
        substr_fasta(fileobj=reference_fasta_obj, chrom=bed_record.chrom,
                     start=ref_start, end=ref_end, out_fa_fn=svp_files_obj.sv_ref_fa)

    if make_subreads_bam:
        # get all raw subreads spanning the structural variants
        make_subreads_bam_of_zmws2(in_subreads_fn_or_obj=subreads_ds_obj, zmws=zmws,
                                   out_bam_fn=svp_files_obj.subreads_bam, out_fa_fn=svp_files_obj.subreads_fa)

    if make_scripts:
        svp_files_obj.make_all_scripts()

    if execute_scripts:
        svp_files_obj.execute_all_scripts(use_sge=use_sge)

    svp_files_obj.make_readme()


def run(args):
    run_polish(in_dir=args.in_dir, in_bed_fn=args.in_bed_fn, out_dir=args.out_dir,
               min_qv=args.min_qv, ref_ext_len=args.ref_ext_len, use_sge=args.use_sge)


def run_polish(in_dir, in_bed_fn, out_dir, min_qv, ref_ext_len, use_sge):
    aln_fn = op.join(in_dir, "alignments.bam")
    subreads_xml_fn = op.join(in_dir, "subreads.xml")
    genome_fa = op.join(in_dir, "genome.fa")
    coverage_fn = op.join(out_dir, 'coverage.txt')

    reference_fasta_obj = Fastafile(genome_fa)
    ofile_obj = open(coverage_fn, 'w')
    alnfile_obj = SingleFileOpener(aln_fn).alignfile
    subreads_ds_obj = SubreadSet(subreads_xml_fn)
    bedreader_obj = BedReader(in_bed_fn)
    samples = bedreader_obj.samples

    report_fp = open(op.join(out_dir, 'report.log'), 'w')

    def report_f(msg):
        report_fp.write(msg + '\n')

    for bed_record, alns in yield_alns_from_bed_file(alnfile_obj, bedreader_obj=bedreader_obj):
        ofile_obj.write('%s\t%s\n' % (len(get_query_subreads_from_alns(alns)),
                                      bed_record.to_str(samples)))  # write coverage info

        sv_prefix = bed2prefix(bed_record)
        work_dir = realpath(op.join(out_dir, sv_prefix))
        log.info("Processing sv %s" % sv_prefix)
        polish_a_sv(bed_record, alns, work_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa=True, make_subreads_bam=True,
                    make_scripts=True, execute_scripts=False, min_qv=min_qv, ref_ext_len=ref_ext_len, use_sge=use_sge, report_f=report_f)

    reference_fasta_obj.close()
    ofile_obj.close()
    alnfile_obj.close()
    subreads_ds_obj.close()
    bedreader_obj.close()
    report_fp.close()


def get_parser():
    """Set up and return argument parser."""
    parser = argparse.ArgumentParser(polish_desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_dir", type=str, help="Input FASTA or FASTQ filename")
    parser.add_argument("in_bed_fn", type=str, help="Structural variants in BED file")
    parser.add_argument("out_dir", type=str, help="Output Directory")
    parser.add_argument("--min_qv", default=C.MIN_POLISH_QV, type=int,
                        help="Minimum Polished QV to include bases in consensus sequence")
    parser.add_argument("--ref_ext_len", default=C.REFERENCE_EXTENSION, type=int,
                        help="Extend reference sequence by ref_ext_len base pairs to both ends")
    parser.add_argument("--use_sge", default=False, action='store_true',
                        help="If True, use SGE; otherwise, run locally")
    return parser


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))


if __name__ == "__main__":
    main()
