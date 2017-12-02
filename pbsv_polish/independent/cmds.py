#!/usr/bin/env python

from pbsv.independent.utils import _is_fmt, cmds_to_bash, execute, realpath, mv_cmd, autofmt, is_fasta
from pbsv.run import svcall_cmd, ngmlrmap_cmd, sort_index_chain_bam_cmd
from .Constants import Constants as C
from .utils import is_fastq
from ..__init__ import TRANSFORM_ENTRY, TRIM_ENTRY, SVDAGCON_ENTRY


def sort_index_bam_inline_cmd(in_fn, nproc=4, tmp_dir=None):
    """Sort and make index of sam or bam files inline"""
    sort_arg = '-T {tmp_dir}/srt' if tmp_dir is not None else ''
    tmp_fn = in_fn + '.sorted'
    c0 = "samtools sort --threads {nproc} {sort_arg} {in_fn} -o {tmp_fn}".format(
        nproc=nproc, sort_arg=sort_arg, in_fn=in_fn, tmp_fn=tmp_fn)
    c1 = mv_cmd(tmp_fn, in_fn)
    c2 = "samtools index {in_fn}".format(in_fn=in_fn)
    return " ; ".join([c0, c1, c2])


def _fn2fmtarg(fn):
    fnext2fmtarg = {"m0": "-m 0", "m4": "-m 4", "bam": "--bam"}
    return fnext2fmtarg[autofmt(fn, fnext2fmtarg.keys())[1]]


def pbindex_cmd(fn):
    return 'pbindex {fn}'.format(fn=fn)


def blasr_cmd(query_fn, target_fn, out_fn, nproc=8):
    return "blasr {q} {t} {fmt} --out {out_fn} --nproc {nproc} --maxMatch 15 --bestn 10 --hitPolicy randombest".\
        format(q=query_fn, t=target_fn, out_fn=out_fn, fmt=_fn2fmtarg(out_fn), nproc=nproc)


def svdagcon_cmd(subreads_bam, output_prefix, output_seq_id, ref_fa):
    return '{} {} {} {} --ref-fa {}'.format(SVDAGCON_ENTRY, subreads_bam, output_prefix, output_seq_id, ref_fa)


def variant_caller_cmd(align_bam, ref_fa, out_fa, out_fq, nproc):
    return "variantCaller --algorithm=best {align_bam} --verbose -j{nproc} --referenceFilename={ref_fa} -o {out_fa} -o {out_fq}".\
        format(align_bam=align_bam, ref_fa=ref_fa, out_fa=out_fa, out_fq=out_fq, nproc=nproc)


def trim_lq_cmd(in_fq, out_fq, out_fa, min_qv):
    assert is_fastq(in_fq) and is_fastq(out_fq) and is_fasta(out_fa)
    # FASTA: simply remove lower case sequences on both ends
    # FASTQ: remove LQ sequences on both ends
    c0 = '{entry} {in_fq} {out_fq} --min-qv {min_qv}'.format(entry=TRIM_ENTRY,
                                                             in_fq=in_fq, out_fq=out_fq, min_qv=min_qv)
    return c0


def sv_transform_coordinate_cmd(in_sv_fn, o_sv_fn):
    return '{entry} {i} {o}'.format(entry=TRANSFORM_ENTRY, i=in_sv_fn, o=o_sv_fn)


def make_input_json_cmd(bam_fn, json_fn, sample):
    """CMD which creates a json file with content [['path_to_bam', 'sample']], which will later
    be used as pbsv call input.
    """
    from pipes import quote
    c0 = 'echo [[\\\"{}\\\", \\\"{}\\\"]] > {}'.format(quote(bam_fn), quote(sample), quote(json_fn))
    return c0

def chain_cmd(in_bam_fn, out_bam_fn, cfg_fn):
    return 'pbsvutil chain {i} {o} --cfg_fn {cfg_fn}'.format(i=in_bam_fn, o=out_bam_fn, cfg_fn=cfg_fn)


def pbsv_run_and_transform_cmds(reads_fn, ref_fa_fn, cfg_fn, o_bam_fn, o_bed_fn, algorithm='ngmlr'):
    """Using ngmlr or blasr to align reads to ref_fa_fn, then run `pbsv call` to call structural variants."""
    o_prefix = o_bam_fn[0:o_bam_fn.rfind('.bam')]
    nochain_bam_fn = o_prefix + '.nochain.bam'
    if algorithm == 'ngmlr':
        # Call pbsvutil ngmlr to map reads_fn to ref_fa_fn to create sorted indexed bam file.
        c0 = ngmlrmap_cmd(in_bam=reads_fn, ref_fn=ref_fa_fn, out_bam=nochain_bam_fn, cfg=cfg_fn)
    elif algorithm == 'blasr':
        c0 = blasr_cmd(query_fn=reads_fn, target_fn=ref_fa_fn, out_fn=nochain_bam_fn)
    else:
        raise ValueError('Could not use algorithm %s to align reads to reference.' % algorithm)

    c1 = chain_cmd(in_bam_fn=nochain_bam_fn, out_bam_fn=o_bam_fn, cfg_fn=cfg_fn)
    c2 = sort_index_bam_inline_cmd(o_bam_fn)

    json_fn = o_prefix + '.json'  # create input json for `pbsv call`
    c3 = make_input_json_cmd(o_bam_fn, json_fn, C.CONSENSUS_SAMPLE)

    tmp_bed = o_bed_fn + '.use_substr_as_chrom.bed'
    c4 = svcall_cmd(ref_fn=ref_fa_fn, in_bam=json_fn, out_bed=tmp_bed, cfg=cfg_fn)
    c5 = rm_ngmlr_indices_cmd(ref_fn=ref_fa_fn)
    c6 = sv_transform_coordinate_cmd(tmp_bed, o_bed_fn)
    return [c0, c1, c2, c3, c4, c5, c6]


def rm_ngmlr_indices_cmd(ref_fn):
    from pbsv.ngmlrmap import _ngm_1, _ngm_2
    return 'rm -f %s %s' % (_ngm_1(ref_fn), _ngm_2(ref_fn))
