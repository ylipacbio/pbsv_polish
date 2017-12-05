#!/usr/bin/env python

import os.path as op
from pbcore.io import SubreadSet, FastaWriter

from pbsv.independent.utils import realpath
from pbsv.cli import _mkdir
from pbsv.io.VcfIO import BedReader, BedWriter, RichBedReader
from pbsv.io.bamstream import SingleFileOpener
from .independent import Constants as C
from .io import SVPolishFiles
from .utils import (Fastafile, bed2prefix, get_query_subreads_from_alns,
        get_query_zmws_from_alns, get_ref_extension_for_sv,
        make_subreads_bam_of_zmws2, substr_fasta, yield_alns_from_bed_file,
        ActionRecord, get_movie_and_zmw_from_name)

import logging
logging.basicConfig(format='%(asctime) %(message)s')
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()

polish_desc = 'Polish input structural variants'

def get_supporting_zmws(sample_reads_dict):
    """
    Return a list of (non-redundant) zmws from sample_reads_dict, where
    sample_reads_dict --- {sample: [reads]}
    ...doctest:
        >>> get_supporting_zmws({'s1': ['movie/0/0_1', 'movie2/2/10_11'], 's2': ['movie/0/13_14']})
        ['movie2/2', 'movie/0']
    """
    zmws = set()
    for sample in sample_reads_dict.keys():
        for r in sample_reads_dict[sample]:
            movie, zmw = get_movie_and_zmw_from_name(r)
            zmws.add('{}/{}'.format(movie, zmw))
    return list(zmws)


def polish_a_sv(bed_record, alns, svobj_dir, subreads_ds_obj, reference_fasta_obj, make_reference_fa, make_subreads_bam, make_scripts, execute_scripts, min_qv, ref_ext_len, use_sge):
    """
    Given a structural variant (as bed_record) and its supportive alignments, polish the given structural variant.
    * if make_reference_fa is True, generate a substr of chromosome as reference for this structural variant
    * if make_subreads_bam is True, generate a subreads bam of all subreads of zmws spanning this strucrtural variant
    * if make_scripts is True, generate scripts to call polished structural variant
    """
    # srs = get_query_subreads_from_alns(alns)
    # zmws = get_query_zmws_from_alns(alns)
    zmws = get_supporting_zmws(bed_record.supporting_reads)
    svp_files_obj = SVPolishFiles(root_dir=svobj_dir, min_qv=min_qv, ref_ext_len=ref_ext_len)
    _mkdir(svp_files_obj.root_dir)  # make a subdirectory (e.g., chrI_0_100_Deletion_-100) for all polishing files
    if make_reference_fa:
        # make a substring spanning the expected structural variants
        ref_start, ref_end = get_ref_extension_for_sv(bed_record,
            reference_fasta_obj.get_reference_length(bed_record.chrom),
            ref_ext_len=ref_ext_len)
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
    return cov >= min_coverage


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
    bedreader_obj = RichBedReader(in_bed_fn)
    bedreader_obj.samples

    skipped_bed_fn, passed_bed_fn = op.join(out_dir, 'in.skipped.bed'), op.join(out_dir, 'in.passed.bed')
    s_writer, p_writer = BedWriter(skipped_bed_fn, bedreader_obj.samples), BedWriter(
        passed_bed_fn, bedreader_obj.samples)

    for bed_record, alns in yield_alns_from_bed_file(alnfile_obj, bedreader_obj=bedreader_obj):
        sv_prefix = bed2prefix(bed_record)
        svobj_dir = realpath(op.join(out_dir, sv_prefix))
        if not sufficient_coverage(bed_record, min_coverage=min_coverage):
            log.info(str(ActionRecord.from_obj(bed_record, action=ActionRecord.PASS, comment='Insufficient_Coverage')))
            s_writer.writeRecord(bed_record)
        else:
            log.info(str(ActionRecord.from_obj(bed_record, action=ActionRecord.PASS, comment='Sufficient_Coverage')))
            p_writer.writeRecord(bed_record)
            polish_a_sv(bed_record, alns, svobj_dir, subreads_ds_obj, reference_fasta_obj,
                        make_reference_fa=True, make_subreads_bam=True,
                        make_scripts=True, execute_scripts=False, min_qv=min_qv, ref_ext_len=ref_ext_len, use_sge=use_sge)

    reference_fasta_obj.close()
    alnfile_obj.close()
    subreads_ds_obj.close()
    bedreader_obj.close()
    s_writer.close()
    p_writer.close()
