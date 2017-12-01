#!/usr/bin/env python

from .independent.Constants import Constants as C
from pbsv.__utils import (get_default_argparser, setup_log, main_runner,
                          compose, subparser_builder, validate_file, args_executer)
from pbsv.run import _mkdir


def add_polish_parser_options(p):
    """Add `pbsvp polish` parser options"""
    fs = [
        _add_in_bed_vcf_parser_option,
        _add_out_dir_parser_option,
        _add_out_bed_vcf_parser_option,
        _add_min_qv_parser_option,
        _add_ref_ext_len_parser_option
    ]
    #f = compose(*fs)
    f = compose(*fs[::-1])
    return f(p)


def add_collect_parser_options(p):
    fs = [
        _add_in_bed_parser_option,
        _add_out_dir_parser_option,
        _add_min_qv_parser_option,
        _add_ref_ext_len_parser_option
    ]
    f = compose(*fs[::-1])


def add_trim_parser_options(p):
    fs = [
        _add_in_fafq_parser_option,
        _add_out_fafq_parser_option,
        _add_qv_windowsize_parser_option,
        _add_min_qv_parser_option
    ]
    f = compose(*fs[::-1])


def add_transform_parser_options(p):
    fs = [
        _add_in_bed_vcf_parser_option,
        _add_out_bed_vcf_parser_option,
    ]
    f = compose(*fs[::-1])


def add_svdagcon_parser_options(p):
    fs = [
        _add_in_subreads_parser_option,
        _add_out_prefix_parser_option,
        _add_consensus_id_parser_option,
        _add_svdagcon_optional_parser_option,
    ]
    f = compose(*fs[::-1])


def _add_in_subreads_parser_option(p):
    p.add_arguments('input_subreads', type=validate_file, help='Input Subreads BAM or SubreadSet')
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


def _add_out_dir_parser_option(p):
    p.add_argument('out_dir', type=_mkdir, help="Output Directory")
    return p


def _add_in_dir_parser_option(p):
    p.add_argument('in_dir', type=str, help="Input Directory")
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


def _add_out_bed_vcf_parser_option(p):
    p.add_argument('out_bed_or_vcf_fn', type=str, help='Output structural variants in BED or VCF')
    return p
