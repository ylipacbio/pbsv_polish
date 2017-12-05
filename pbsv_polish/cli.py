#!/usr/bin/env python

"""
Define command line entry for tool `pbsvp`, including
    pbsvp polish
    pbsvp collect

Define `pbsvputil`, utils for `pbsvp`, including
    pbsvputil trim-lq
    pbsvputil transform-coordinate
    pbsvputil svdagcon
"""

import sys
import logging
from .__init__ import (get_version, POLISH_ENTRY, COLLECT_ENTRY, TRIM_ENTRY,
                       SVDAGCON_ENTRY, TRANSFORM_ENTRY)
from pbsv.__utils import (get_default_argparser, setup_log, main_runner,
                          compose, subparser_builder, validate_file, args_executer)
from .argsutil import (add_polish_parser_options,  add_collect_parser_options,
                       add_trim_parser_options, add_transform_parser_options, add_svdagcon_parser_options)
from .polish import polish_desc, run_polish
from .collect import collect_desc, run_collect
from .trim_lq import trim_desc, run_trim
from .transform_coordinates import transform_desc, run_transform
from .svdagcon import svdagcon_desc, run_svdagcon


log = logging.getLogger()
slog = logging.getLogger('status.' + __file__)


def _args_run_polish(args):
    """Run `pbsvp polish`"""
    log.info("Running `{}`".format(POLISH_ENTRY))
    log.debug('Locals={}'.format(locals()))
    run_polish(genome_fa=args.genome_fa, subreads_xml_fn=args.subreads_bam, aln_fn=args.alignments_bam,
               in_bed_fn=args.in_rich_bed, out_dir=args.out_dir,
               min_coverage=args.min_coverage, min_qv=args.min_qv,
               ref_ext_len=args.ref_ext_len, use_sge=args.use_sge)
    return 0


def _args_run_collect(args):
    """Run `pbsvp collect`"""
    log.info("Running `{}`".format(COLLECT_ENTRY))
    log.debug('Locals={}'.format(locals()))
    run_collect(work_dir=args.work_dir, collected_bed_fn=args.out_bed_or_vcf_fn,
                min_qv=args.min_qv, ref_ext_len=args.ref_ext_len)
    return 0


def _args_run_trim(args):
    log.info("Running `{}`".format(TRIM_ENTRY))
    log.debug('Locals={}'.format(locals()))
    run_trim(i_fn=args.in_fa_or_fq_fn, o_fn=args.out_fa_or_fq_fn,
             min_qv=args.min_qv, windowsize=args.qv_windowsize)
    return 0


def _args_run_svdagcon(args):
    log.info("Running `{}`".format(SVDAGCON_ENTRY))
    log.debug('Locals={}'.format(locals()))
    run_svdagcon(input_subreads_bam=args.subreads_bam, ref_fa=args.ref_fa,
                 output_prefix=args.output_prefix, consensus_id=args.consensus_id,
                 nproc=args.nproc, max_score=args.max_score,
                 use_first_seq_if_fail=args.use_first_seq_if_fail)
    return 0


def _args_run_transform(args):
    log.info('Running `{}`'.format(TRANSFORM_ENTRY))
    log.debug('Locals={}'.format(locals()))
    run_transform(i_fn=args.in_bed_or_vcf_fn, o_fn=args.out_bed_or_vcf_fn)
    return 0


def pbsvp_get_parser():
    """Get parser for pbsvp subcommands"""
    desc = "PacBio Structural Variants Polish Tool Suite"
    p = get_default_argparser(version=get_version(), description=desc)
    sp = p.add_subparsers(help='commands')

    def builder(subparser_id, description, options_func, exe_func, epilog=None):
        """subparser builder"""
        subparser_builder(sp, subparser_id, description, options_func, exe_func, epilog)

    # `pbsvp polish`, polish structural variants
    polish_desc = "Polish structural variants"
    builder('polish', polish_desc, add_polish_parser_options, _args_run_polish)

    # `pbsvp collect`, collect polished structural variants
    collect_desc = "Collect polished structural variants"
    builder('collect', collect_desc, add_collect_parser_options, _args_run_collect)
    return p


def pbsvputil_get_parser():
    """Get parser for pbsvputil subcommands"""
    desc = "PacBio Structural Variants Polish Utils"
    p = get_default_argparser(version=get_version(), description=desc)
    sp = p.add_subparsers(help='commands')

    def builder(subparser_id, description, options_func, exe_func, epilog=None):
        """subparser builder"""
        subparser_builder(sp, subparser_id, description, options_func, exe_func, epilog)

    # `pbsvputil trim-lq`
    builder('trim-lq', trim_desc, add_trim_parser_options, _args_run_trim)

    # `pbspvutil svdagcon`
    builder('svdagcon', svdagcon_desc, add_svdagcon_parser_options, _args_run_svdagcon)

    # `pbsvputil transform-coordinate`
    builder('transform-coordinate', transform_desc, add_transform_parser_options, _args_run_transform)
    return p


def pbsvp_main(argv=None):
    """pbsvp Main function, entry for command line tool `pbsvp`"""
    argv_ = sys.argv if argv is None else argv
    parser = pbsvp_get_parser()
    return main_runner(argv_[1:], parser, args_executer, setup_log, log)


def pbsvputil_main(argv=None):
    """pbsvputil main function, entry for command line tool `pbsvputil`"""
    argv_ = sys.argv if argv is None else argv
    parser = pbsvputil_get_parser()
    return main_runner(argv_[1:], parser, args_executer, setup_log, log)
