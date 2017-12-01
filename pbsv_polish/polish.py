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
from .utils import BedReader, Fastafile, SubreadSet, bed2prefix, get_query_subreads_from_alns, get_query_zmws_from_alns, get_ref_extension_for_sv, make_subreads_bam_of_zmws2, substr_fasta, yield_alns_from_bed_file

import logging
logging.basicConfig(format='%(asctime) %(message)s')
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()

polish_desc = 'Polish input structural variants'


def polish_a_sv(bed_record, alns, work_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa, make_subreads_bam, make_scripts, execute_scripts, min_qv, ref_ext_len, use_sge):
    """
    Given a structural variant (as bed_record) and its supportive alignments, polish the given structural variant.
    * if make_reference_fa is True, generate a substr of chromosome as reference for this structural variant
    * if make_subreads_bam is True, generate a subreads bam of all subreads of zmws spanning this strucrtural variant
    * if make_scripts is True, generate scripts to call polished structural variant
    """
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


def sufficient_coverage(bed_record, min_coverage):
    """Return True if bed_record has sufficient coverage of supportive reads,
    otherwise, False"""
    cov = sum([bed_record.fmts[sample].ad for sample in bed_record.samples])
    return cov <= min_coverage


class ActionRecord(object):
    SEP = '\t'
    HEADER = SEP.join(['name', 'action', 'ad', 'dp', 'comment'])
    SKIP, PASS, REJECTED, POLIHSED = ['SKIP', 'PASS', 'REJECTED', 'POLISHED']

    def __init__(self, name, action, ad, dp, comment):
        self.name, self.action, self.ad, self.dp, self.comment = name, action, int(ad), int(dp), comment

    def __str__(self):
        return self.SEP.join([str(s) for s in [self.name, self.action, self.ad, self.dp, self.comment]])

    @classmethod
    def from_obj(cls, svobj, action, comment):
        ad = sum([svobj.fmts[sample].ad for sample in svobj.samples])
        dp = sum([svobj.fmts[sample].dp for sample in svobj.samples])
        return ActionRecord(bed2prefix(svobj), action, ad, dp, comment)

    @classmethod
    def from_str(cls, line):
        name, action, ad, dp, comment = line.split(cls.SEP)[0:5]
        return ActionRecord(name, action, ad, dp, comment)


def run_polish(genome_fa, subreads_xml_fn, aln_fn, in_bed_fn, out_dir,
               min_coverage, min_qv, ref_ext_len, use_sge):
    """
    Given input genome FASTA (genome_fa), PACBIO subreads (subreads_xml or bam),
    alignments bam mapping subreads to genome, pbsv output structural variants
    (in_bed_fn), polish structural variants.
    min_coverage --- Skip structural variants not having enough coverage of
         supportive reads
    ref_ext_len --- Extend corresponding reference substrings of structural variants
         to both ends by ref_ext_len base pairs.
    min_qv --- Trim low quality endings of structural variants consensus sequence.
    """
    reference_fasta_obj = Fastafile(genome_fa)
    alnfile_obj = SingleFileOpener(aln_fn).alignfile
    subreads_ds_obj = SubreadSet(subreads_xml_fn)
    bedreader_obj = BedReader(in_bed_fn)
    samples = bedreader_obj.samples

    action_record_fp = open(op.join(out_dir, 'actions.log'), 'w')
    action_record_fp.write(ActionRecord.HEADER + '\n')

    for bed_record, alns in yield_alns_from_bed_file(alnfile_obj, bedreader_obj=bedreader_obj):
        sv_prefix = bed2prefix(bed_record)
        work_dir = realpath(op.join(out_dir, sv_prefix))
        log.info("Processing sv %s" % sv_prefix)

        if not sufficient_coverage(bed_record, min_coverage=min_coverage):
            ar = ActionRecord.from_obj(bed_record, action=ActionRecord.SKIP, comment='Insufficient Coverage')
            continue
        else:
            ar = ActionRecord.from_obj(bed_record, action=ActionRecord.PASS, comment='Sufficient Coverage')
        action_record_fp.write(str(ar) + '\n')
        polish_a_sv(bed_record, alns, work_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa=True, make_subreads_bam=True,
                    make_scripts=True, execute_scripts=False, min_qv=min_qv, ref_ext_len=ref_ext_len, use_sge=use_sge)

    reference_fasta_obj.close()
    alnfile_obj.close()
    subreads_ds_obj.close()
    bedreader_obj.close()
    action_record_fp.close()
