#!/usr/bin/env python
import os
import os.path as op
from pbsv.__utils import (get_default_argparser, setup_log, main_runner,
                          compose, subparser_builder, args_executer)
from .independent.Constants import Constants as C
from .independent.utils import get_zmws_from_iter


def mkdir(path):
    if not op.isdir(path):
        os.makedirs(path)
    return path


def validate_file(fn):
    if op.exists(fn):
        return fn

def parse_zmws_from_str(s, sep=';'):
    """
    ...doctest:
        >>> parse_zmws_from_str('movie/1;movie/2;movie2/3')
        ['movie/1', 'movie/2', 'movie2/3']
    """
    return get_zmws_from_iter(s.split(sep))


def add_polish_parser_options(p):
    """Add `pbsvp polish` parser options"""
    fs = [
        _add_genome_fa_parser_option,
        _add_subreads_bam_parser_option,
        _add_alignment_bam_parser_option,
        _add_in_rich_bed_parser_option,
        #_add_in_bed_vcf_parser_option,
        _add_out_dir_parser_option,
        _add_min_coverage_parser_option,
        _add_min_qv_parser_option,
        _add_ref_ext_len_parser_option,
        _add_use_sge_parser_option
    ]
    f = compose(*fs[::-1])
    return f(p)


def add_collect_parser_options(p):
    fs = [
        _add_work_dir_parser_option,
        _add_out_bed_vcf_parser_option,
        _add_min_qv_parser_option,
        _add_ref_ext_len_parser_option
    ]
    f = compose(*fs[::-1])
    return f(p)


def add_trim_parser_options(p):
    fs = [
        _add_in_fafq_parser_option,
        _add_out_fafq_parser_option,
        _add_qv_windowsize_parser_option,
        _add_min_qv_parser_option
    ]
    f = compose(*fs[::-1])
    return f(p)


def add_transform_parser_options(p):
    fs = [
        _add_in_bed_vcf_parser_option,
        _add_out_bed_vcf_parser_option,
    ]
    f = compose(*fs[::-1])
    return f(p)


def add_svdagcon_parser_options(p):
    fs = [
        _add_subreads_bam_parser_option,
        _add_out_prefix_parser_option,
        _add_consensus_id_parser_option,
        _add_svdagcon_optional_parser_option
    ]
    f = compose(*fs[::-1])
    return f(p)


def add_extractsr_parser_options(p):
    fs = [
        _add_in_subreads_bam_parser_option,
        _add_zmws_parser_option,
        _add_out_subreads_bam_parser_option,
    ]
    f = compose(*fs[::-1])
    return f(p)

def _add_zmws_parser_option(p):
    p.add_argument('zmws', type=parse_zmws_from_str, help='semicolon separated zmws, e.g., movie/zmw1;movie/zmw2;movie2/zmw3')
    return p


def _add_in_subreads_bam_parser_option(p):
    p.add_argument('in_subreads_bam', type=validate_file, help='Input Subreads BAM or SubreadSet')
    return p


def _add_out_subreads_bam_parser_option(p):
    p.add_argument('out_subreads_bam', type=validate_file, help='Output Subreads BAM or SubreadSet')
    return p


def _add_subreads_bam_parser_option(p):
    p.add_argument('subreads_bam', type=validate_file, help='Subreads BAM or SubreadSet')
    return p


def _add_alignment_bam_parser_option(p):
    p.add_argument('alignments_bam', type=validate_file, help='Alignments BAM')
    return p


def _add_out_prefix_parser_option(p):
    p.add_argument("output_prefix", help="Output filename prefix (ex: sv_consensus)")
    return p


def _add_consensus_id_parser_option(p):
    p.add_argument("consensus_id", help="Consensus sequence ID name (ex: chr1_100_100_Insertion)")
    return p


def _add_svdagcon_optional_parser_option(p):
    p.add_argument("--nproc", default=8, type=int, help="Number of processes")
    p.add_argument("--max-score", default=-1000, type=int, help="Maximum score to use a blasr alignment")
    p.add_argument("--use-first-seq-if-fail", default=True, action='store_false',
                   help="Use the first sequence as backup reference if pbdagcon fails")
    p.add_argument("--ref-fa", type=str,
                   help="Use reference fasta to bound consensus sequence and remove unmappablei edges")
    return p


def _add_in_bed_parser_option(p):
    p.add_argument('in_bed_fn', type=validate_file, help="Input structural variants in BED")
    return p


def _add_work_dir_parser_option(p):
    p.add_argument('work_dir', type=str, help="Working Directory")
    return p


def _add_out_dir_parser_option(p):
    p.add_argument('out_dir', type=mkdir, help="Output Directory")
    return p


def _add_genome_fa_parser_option(p):
    p.add_argument('genome_fa', type=validate_file, help="Reference Genome FASTA")
    return p


def _add_min_coverage_parser_option(p):
    p.add_argument('--min-coverage', type=int, default=C.MIN_POLISH_COVERAGE,
                   help="Minimum number of supportive reads to polish a strucutural variant")
    return p


def _add_min_qv_parser_option(p):
    p.add_argument("--min-qv", default=C.MIN_POLISH_QV, type=int,
                   help="Minimum Polished QV to include bases in consensus sequence")
    return p


def _add_ref_ext_len_parser_option(p):
    p.add_argument("--ref-ext-len", default=C.REFERENCE_EXTENSION, type=int,
                   help="Extend reference sequence by ref_ext_len base pairs to both ends")
    return p


def _add_use_sge_parser_option(p):
    p.add_argument("--use-sge", default=False, action='store_true',
                   help="If True, use SGE; otherwise, run locally")
    return p


def _add_qv_windowsize_parser_option(p):
    p.add_argument('--qv-windowsize', default=100, type=int,
                   help='Compute average MapQV in windows of size')
    return p


def _add_in_fafq_parser_option(p):
    p.add_argument('in_fa_or_fq_fn', type=validate_file, help='Input FASTA or FASTQ file')
    return p


def _add_out_fafq_parser_option(p):
    p.add_argument('out_fa_or_fq_fn', type=str, help='Output FASTA or FASTQ file')
    return p


def _add_in_bed_vcf_parser_option(p):
    p.add_argument('in_bed_or_vcf_fn', type=validate_file, help='Input structural variants in BED or VCF')
    return p

def _add_in_rich_bed_parser_option(p):
    p.add_argument('in_rich_bed', type=validate_file, help='Input structural variants in rich BED.')
    return p

def _add_out_bed_vcf_parser_option(p):
    p.add_argument('out_bed_or_vcf_fn', type=str, help='Output structural variants in BED or VCF')
    return p
