#!/usr/bin/env python

from argparse import ArgumentParser
import sys

from .independent.utils import is_fasta, is_fastq
from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter

trim_desc = "Trim LQ sequences on both ends, where average MapQV of LQ sequences are less than min_qv." + \
            "If Output is out.fastq, also write to out.fasta."

def get_hq_start_end(seq, is_hq_f):
    """
    ..doctest:
    >>> def f(c): return c.isupper()
    >>> get_hq_start_end('', f)
    (0, 0)
    >>> get_hq_start_end('GC', f)
    (0, 2)
    >>> get_hq_start_end('aGCc', f)
    (1, 3)
    """
    start, end = 0, 0
    i = 0
    while i < len(seq):
        if is_hq_f(seq[i]):
            start = i
            break
        i += 1

    i = len(seq) - 1
    while i >= 0:
        if is_hq_f(seq[i]):
            end = i + 1
            break
        i -= 1
    return (start, end)


def get_hq_start_end_fasta(seq):
    def f(c):
        return c.isupper()
    return get_hq_start_end(seq, f)


def get_hq_start_end_fastq(qual, min_qv):
    def f(c):
        return c >= min_qv
    return get_hq_start_end(qual, f)


def trim_fasta(i_fn, o_fn):
    with FastaReader(i_fn) as reader, FastaWriter(o_fn) as writer:
        for r in reader:
            hq_start, hq_end = get_hq_start_end_fasta(r.sequence)
            writer.writeRecord(r.name + '/%s_%s' % (hq_start, hq_end), r.sequence[hq_start:hq_end])


def get_avg_value_in_window(values, windowsize):
    """Given a list of integers, return average value in windows, where ret[i] is the average value of
    values[i, min(i+windowsize, len(values))]
    """
    ret = [0] * len(values)
    for i in range(0, len(values)):
        end = min(i + windowsize, len(values))
        ret[i] = sum(values[i:end]) * 1.0 / (end - i)
    return ret


def trim_fastq(i_fn, o_fn, o_fa_fn, min_qv, windowsize):
    with FastqReader(i_fn) as reader, FastqWriter(o_fn) as writer, FastaWriter(o_fa_fn) as fawriter:
        for r in reader:
            avg_qual = get_avg_value_in_window(r.quality, windowsize)
            hq_start, hq_end = get_hq_start_end_fastq(avg_qual, min_qv)
            newid = r.id + '/%s_%s' % (hq_start, hq_end)
            newseq = r.sequence[hq_start:hq_end]
            newqul = r.quality[hq_start:hq_end]
            writer.writeRecord(newid, newseq, newqul)
            fawriter.writeRecord(newid, newseq)


def run_trim(i_fn, o_fn, min_qv, windowsize):
    if all(is_fasta(fn) for fn in [i_fn, o_fn]):
        trim_fasta(i_fn, o_fn)
    elif all(is_fastq(fn) for fn in [i_fn, o_fn]):
        o_fa_fn = o_fn[0:o_fn.rfind('.')] + '.fasta'
        trim_fastq(i_fn, o_fn, o_fa_fn, min_qv, windowsize)
    else:
        raise ValueError("Input and output must be both BED or VCF")
