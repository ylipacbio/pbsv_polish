#!/usr/bin/env python

"""
Make a consensus sequence from subreads.bam
(1) Find the best seed as reference
(2) Align rest to seed
(3) Call pbdagcon
"""

from argparse import ArgumentParser
from collections import defaultdict
import os
import os.path as op
import sys
import json
import numpy as np

from pbcore.util.Process import backticks
from pbcore.io import DataSet, FastaReader, FastaWriter
from pbsv.ngmlrmap import ln_cmd
from pbsv.independent.utils import execute, realpath
from pbsv.cli import _mkdir
from pbsv.io.linefile import X2PysamReader
from pbsv.libs import AlignmentFile, AlignedSegment

from pbtranscript.io import FastaRandomReader, BLASRM4Reader
from pbsv.__init__ import get_version


class AlignGraphUtilError(Exception):
    """Align Group Util Error Class"""
    pass


def choose_template_by_blasr(fasta_filename, out_filename, nproc=8,
                             maxScore=-1000, min_number_reads=1):
    """
    Choose the best template for gcon reference
    Pick the one that has the highest sum of score to others

    Returns: FastaRecord of selected ref
    """
    fd = FastaRandomReader(fasta_filename)

    cmd = "blasr --nproc {nproc} ".format(nproc=nproc) + \
          "--maxScore {score} ".format(score=maxScore) + \
          "--maxLCPLength 15 --bestn 10 --nCandidates 50 " + \
          "-m 1 {fa} {fa} ".format(fa=fasta_filename) + \
          "--out {out} ".format(out=out_filename) + \
          "1>/dev/null 2>/dev/null"

    out, code, msg = backticks(cmd)
    if code != 0:
        return None

    # blasr -m 1 output format:
    # (0) qName   (1) tName (2) qStrand
    # (3) tStrand (4) score (5) percentSimilarity
    # (6) tStart  (7) tEnd  (8) tLength
    # (9) qStart  (10) qEnd (11) qLength
    # (12) nCells
    scores = defaultdict(lambda: [])
    with open(out_filename) as f:
        for line in f:
            raw = line.strip().split()
            # qID gets an extra /0_length
            qID, tID = raw[0][:raw[0].rfind('/')], raw[1]
            if qID == tID:
                continue  # self-hit, ignore
            scores[qID].append(float(raw[4]))  # use score

    # find the one with the highest average alignment similarity
    #score_array = []
    #for k, v in scores.iteritems():
    #     score_array.append((np.ceil(np.mean(v)), k))

    #if len(score_array) < min_number_reads:
    if len(scores) < min_number_reads:
        errMsg = "Not enough number of reads in " + \
                 "choose_template_by_blasr {0} < {1}".format(
                     len(scores), min_number_reads)
        raise AlignGraphUtilError(errMsg)

    #score_array.sort(reverse=True)

    # Find the longest sequence that is within the std deviation of
    #best_mean, best_id = score_array[0]
    #best_len = len(fd[best_id].sequence)
    #for _mean, _id in score_array[1:]:
    #    if _mean != best_mean:
    #        break
    #    _len = len(fd[_id].sequence)
    #    if _len > best_len:
    #        best_id = _id
    #        best_len = _len
    best_id, best_score = None, 0
    for k, v in scores.iteritems():
        this_score = sum(v)
        if this_score < best_score: # blasr score the less the better
            best_id = k
            best_score = this_score

    return fd[best_id]


def make_aln_input_to_ref(fasta_filename, ref_filename,
                          out_filename, nproc=8):
    """
    Make blasr -m 5 output of input aligned to ref

    Return: alignment filename
    """
    # pbdagcon only takes -m 5 output
    tmp_out = "{out}.tmp".format(out=out_filename)

    cmd = "blasr {infa} ".format(infa=fasta_filename) + \
          "{ref} --bestn 1 ".format(ref=ref_filename) + \
          "--nproc {nproc} ".format(nproc=nproc) + \
          "-m 5 --out {out} ".format(out=tmp_out) + \
          "1>/dev/null 2>/dev/null"

    out, code, msg = backticks(cmd)
    if code != 0:
        errMsg = "Unable to align {infa} to {ref} in make_aln_input_to_ref".\
            format(infa=fasta_filename, ref=ref_filename)
        raise AlignGraphUtilError(errMsg)

    # trim away anything that is a self-hit or opp strand
    with open(out_filename, 'w') as f, \
            open(tmp_out, 'r') as h:
        for line in h:
            raw = line.strip().split()
            # blasr -m 5 output format:
            # (0) qName (1) qLength (2) qStart (3) qEnd (4) qStrand
            # (5) tName (6) tLength (7) tStart (8) tEnd (9) tStrand
            # (10...) score ...
            if raw[0] == raw[5]:
                continue  # self-hit
            if raw[4] != raw[9]:
                continue  # opp strand
            f.write(line)

    os.remove(tmp_out)


def pbdagcon_wrapper(fasta_filename, output_prefix, consensus_name, nproc=8,
                     maxScore=-1000, min_seq_len=300, use_first_seq_if_fail=True):
    """
    (1) Find the best seed as reference
    (2) Align rest to seed
    (3) Call pbdagcon
    """
    ref_filename = output_prefix + '_ref.fasta'
    cons_filename = output_prefix + '.fasta'
    tmp_cons_filename = output_prefix + '.fasta.tmp'
    aln_filename = output_prefix + '.saln'
    try:
        out_filename_m1 = output_prefix + ".saln.m1"
        ref = choose_template_by_blasr(fasta_filename=fasta_filename,
                                       out_filename=out_filename_m1,
                                       nproc=nproc, maxScore=maxScore)
        os.remove(out_filename_m1)

        with open(ref_filename, 'w') as f:
            f.write(">{0}\n{1}\n".format(consensus_name, ref.sequence))

        # create alignment file
        make_aln_input_to_ref(fasta_filename=fasta_filename,
                              ref_filename=ref_filename,
                              out_filename=aln_filename,
                              nproc=nproc)
        # call pbdagcon
        cmd = "pbdagcon -t 0 -m {minlen} -c 1 -j {nproc} {aln} > {out}".format(
            minlen=min_seq_len, nproc=nproc, aln=aln_filename,
            out=tmp_cons_filename)
        cmd += " 2>/dev/null"
        out, code, msg = backticks(cmd)
        if code != 0:
            raise AlignGraphUtilError("Cannot run command: %s" % cmd)

        with FastaReader(tmp_cons_filename) as reader, \
            open(cons_filename, 'w') as writer:
            for rec in reader:
                name = rec.name.strip()
                if "/" in name:
                    # change cid format from c{cid}/0_{len} to c{cid}
                    name = name[:name.find('/')]
                seq = rec.sequence.strip()
                if not 'N' in seq: # Don't write if seq contains N
                    writer.write(">{0}\n{1}\n".format(name, seq))
        os.remove(tmp_cons_filename)
    except AlignGraphUtilError:
        # pick the first sequence as reference as a backup plan
        if use_first_seq_if_fail:
            first_seq = FastaReader(fasta_filename).__iter__().next()
            with open(ref_filename, 'w') as f:
                f.write(">{0}_ref\n{1}\n".
                        format(consensus_name, first_seq.sequence))
        else:
            #("Could not make pbdgacon consensus for %s" % fasta_filename)
            with open(ref_filename, 'w') as f:
                f.write("") # empty output
        if op.exists(tmp_cons_filename):
            os.remove(tmp_cons_filename)
    return 0

def make_fai(fn):
    execute('samtools faidx %s' % fn)

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser()
    parser.add_argument("input_subreads_bam", help="Input fasta filename")
    parser.add_argument("output_prefix", help="Output filename prefix (ex: sv_consensus)")
    parser.add_argument("consensus_id", help="Consensus sequence ID name (ex: chr1_100_100_Insertion)")
    parser.add_argument("--nproc", default=8, type=int, help="Number of processes")
    parser.add_argument("--maxScore", default=-1000, type=int, help="blasr maxScore")
    parser.add_argument("--use_first_seq_if_fail", default=True, action='store_false', help="Use the first sequence as backup reference if pbdagcon fails")
    parser.add_argument("--ref_fa", type=str, help="Use reference fasta to bound consensus sequence and remove unmappablei edges")
    parser.add_argument("--version", "-v", action='version', version='%(prog)s ' + get_version())
    return parser


def restore_args_with_whitespace(x):
    """Restore args with whitespace, to support file paths with whitespace"""
    y = []
    for xi in x:
        if len(y) == 0:
            y.append(xi)
        else:
            if y[-1].endswith('\\'):
                y[-1] = y[-1][0:-1] + ' ' + xi
            else:
                y.append(xi)
    return y

def get_fasta_fn_from_subreads_bam_fn(bam_fn):
    fasta_fn = op.join(op.dirname(bam_fn), op.basename(bam_fn).split('.')[0] + '.subreads.fasta')
    alnfile = X2PysamReader(bam_fn)._alignment_file
    writer = FastaWriter(fasta_fn)
    for aln in alnfile:
        writer.writeRecord(aln.query_name, aln.query_sequence)
    writer.close()
    return fasta_fn

from .utils import blasr_cmd

def write_best_mappable_substr(in_fa_fn, ref_fa_fn, out_fa_fn, read_name):
    from .utils import substr_fasta
    from pbsv.libs import Fastafile

    forward_m4_fn = out_fa_fn+'.forward.m4'
    reverse_m4_fn = out_fa_fn+'.reverse.m4'
    execute(blasr_cmd(query_fn=in_fa_fn, target_fn=ref_fa_fn, out_fn=forward_m4_fn))
    execute(blasr_cmd(query_fn=ref_fa_fn, target_fn=in_fa_fn, out_fn=reverse_m4_fn))
    start, end = get_span_region_from_m4(forward_m4_fn, reverse_m4_fn)
    if start is not None and end is not None:
        substr_fasta(fileobj=Fastafile(in_fa_fn), chrom=read_name, start=start, end=end, o_fasta_fn=out_fa_fn)
    else: # otherwise, output is empty
        open(out_fa_fn, 'w').write('')


def get_span_region_from_m4(forward_m4_fn, reverse_m4_fn):
    fs, fe = _get_query_span_region_from_m4(forward_m4_fn)
    rs, re = _get_target_span_region_from_m4(reverse_m4_fn)
    def f(a, b, op):
        if (a is None and b is None):
            return None
        elif a is None:
            return b
        elif b is None:
            return a
        else:
            return op(a, b)
    return f(fs, rs, min), f(fe, re, max)

def _get_query_span_region_from_m4(m4_fn):
    start, end = None, None
    for r in BLASRM4Reader(m4_fn):
        start = r.qStart if start is None else min(start, r.qStart)
        end = r.qEnd if end is None else max(end, r.qEnd)
    return start, end

def _get_target_span_region_from_m4(m4_fn):
    start, end = None, None
    for r in BLASRM4Reader(m4_fn):
        if r.sStrand == '+':
            start = r.sStart if start is None else min(start, r.sStart)
            end = r.sEnd if end is None else max(end, r.sEnd)
        elif r.sStrand == '-':
            start = (r.sLength-r.sEnd) if start is None else min(start, (r.sLength-r.sEnd)) # exclusive end to inclusive start
            end = (r.sLength-r.sStart) if end is None else max(end, (r.sLength-r.sStart)) # inclusive start to exclusive end
    return start, end


def run(args):
    """Build consesus sequences from subreads.bam input using pbdagcon."""
    #convert subreads.bam to subread.fasta
    sr_fasta = get_fasta_fn_from_subreads_bam_fn(args.input_subreads_bam)
    op.exists(sr_fasta)
    pbdagcon_wrapper(fasta_filename=sr_fasta, output_prefix=args.output_prefix+'.all', consensus_name=args.consensus_id, nproc=args.nproc, maxScore=args.maxScore, use_first_seq_if_fail=args.use_first_seq_if_fail)
    if args.ref_fa:
        write_best_mappable_substr(in_fa_fn=args.output_prefix+'.all.fasta', ref_fa_fn=args.ref_fa, out_fa_fn=args.output_prefix+'.fasta', read_name=args.consensus_id)
    else:
        execute('ln %s %s' % (args.output_prefix+'.all.fasta', args.output_prefix+'.fasta'))

    try: # OK to fail indexing fasta, this may happen if fasta is empty
        make_fai(args.output_prefix+'.fasta')
    except Exception:
        pass


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
