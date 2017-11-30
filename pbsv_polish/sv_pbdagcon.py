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
import numpy as np

from pbcore.util.Process import backticks
from pbcore.io import FastaReader, FastaWriter
from pbsv.independent.utils import execute, realpath
from pbsv.cli import _mkdir
from pbsv.io.bamstream import SingleFileOpener
from pbsv.libs import AlignmentFile, Fastafile

from pbtranscript.io import FastaRandomReader, BLASRM4Reader
from pbsv.__init__ import get_version
from .utils import *


class AlignGraphUtilError(Exception):
    """Align Group Util Error Class"""


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

    out, code, _ = backticks(cmd)
    if code != 0:
        return None

    # blasr -m 1 output format:
    # (0) qID   (1) tID (2) qStrand
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
            # scores[qID].append(float(raw[4]))  # use score
            scores[qID].append(float(raw[5]) / 100.0 * (int(raw[10]) - int(raw[9]))
                               )  # use similarity * aligned query length

    # find the one with the highest average alignment similarity
    #score_array = []
    # for k, v in scores.iteritems():
    #     score_array.append((np.ceil(np.mean(v)), k))

    # if len(score_array) < min_number_reads:
    if len(scores) < min_number_reads:
        errMsg = "Not enough number of reads in " + \
                 "choose_template_by_blasr {0} < {1}".format(
                     len(scores), min_number_reads)
        raise AlignGraphUtilError(errMsg)

    # score_array.sort(reverse=True)

    # Find the longest sequence that is within the std deviation of
    #best_mean, best_id = score_array[0]
    #best_len = len(fd[best_id].sequence)
    # for _mean, _id in score_array[1:]:
    #    if _mean != best_mean:
    #        break
    #    _len = len(fd[_id].sequence)
    #    if _len > best_len:
    #        best_id = _id
    #        best_len = _len
    def get_sr_per_zmw(fasta_filename):
        """Given a subreads fasta file, return number of subreads in each zmw"""
        sr_per_zmw = defaultdict(lambda: 0)
        from pbsv.libs import Fastafile
        from .utils import zmw_from_subread
        obj = Fastafile(fasta_filename)
        for name in obj.references:
            sr_per_zmw[zmw_from_subread(name)] += 1
        return sr_per_zmw

    sr_per_zmw = get_sr_per_zmw(fasta_filename)
    best_id, best_score = None, 0
    for k, v in scores.iteritems():
        sr_per_zmw[zmw_from_subread(k)]
        this_score = sum(v)
        # if this_score < best_score: # blasr score the less the better
        if this_score > best_score:  # similarity * aligned length, the larger the better
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

    out, code, _ = backticks(cmd)
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
            # (0) qID (1) qLength (2) qStart (3) qEnd (4) qStrand
            # (5) sID (6) tLength (7) tStart (8) tEnd (9) tStrand
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
        out, code, _ = backticks(cmd)
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
                if not 'N' in seq:  # Don't write if seq contains N
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
                f.write("")  # empty output
        if op.exists(tmp_cons_filename):
            os.remove(tmp_cons_filename)
    return 0


def make_fai(fn):
    """Make FAI index file of a fasta file."""
    assert ' ' not in fn
    execute('rm -f %s' % (fn + '.fai'))
    execute('samtools faidx %s' % fn)


def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser()
    parser.add_argument("input_subreads_bam", help="Input fasta filename")
    parser.add_argument("output_prefix", help="Output filename prefix (ex: sv_consensus)")
    parser.add_argument("consensus_id", help="Consensus sequence ID name (ex: chr1_100_100_Insertion)")
    parser.add_argument("--nproc", default=8, type=int, help="Number of processes")
    parser.add_argument("--maxScore", default=-1000, type=int, help="blasr maxScore")
    parser.add_argument("--use_first_seq_if_fail", default=True, action='store_false',
                        help="Use the first sequence as backup reference if pbdagcon fails")
    parser.add_argument("--ref_fa", type=str,
                        help="Use reference fasta to bound consensus sequence and remove unmappablei edges")
    parser.add_argument("--version", "-v", action='version', version='%(prog)s ' + get_version())
    return parser


def get_fasta_fn_from_subreads_bam_fn(bam_fn):
    """
    Given a bam file `*.bam`, write record to fasta in the same directory as `*.subreads.fasta`
    """
    fasta_fn = op.join(op.dirname(bam_fn), op.basename(bam_fn).split('.')[0] + '.subreads.fasta')
    alnfile = SingleFileOpener(fn=bam_fn).alignfile
    writer = FastaWriter(fasta_fn)
    for aln in alnfile:
        writer.writeRecord(aln.query_name, aln.query_sequence)
    writer.close()
    return fasta_fn


from .utils import blasr_cmd


def get_region_of_seq_in_a_match_b(a_fa_obj, b_fa_obj, a_seq_name, work_dir):
    """
    a_fa_obj --- Fastafile a, which must contain sequence a_seq_name, and may contain other sequences
    b_fa_obj --- Fastafile b, which contains exactly one sequence
    a_seq_name -- name of a sequence which a_fa_obj must contain
    work_dir --- working directory for all intermediate files.

    Return (chrom, start, end) of a region which
        * is a substring of a_seq_name sequence of a_fa_obj, and
        * covers alignments that match the only sequence in b_fa_obj
    and return (a_seq_name, start, end).
    """
    assert isinstance(a_fa_obj, Fastafile) and isinstance(b_fa_obj, Fastafile)
    a_fa_fn, b_fa_fn = a_fa_obj.filename, b_fa_obj.filename
    if a_seq_name and not (a_seq_name in a_fa_obj.references):
        raise ValueError("FASTA file %s does not contain sequence %s" % (a_fa_fn, a_seq_name))
    assert_fasta_has_one_seq(b_fa_obj)
    # call `blasr` to align a and b forward and backward (aka a to b and b to a).
    a2b_m4_fn = op.join(work_dir, '%s.vs.%s.m4' % (basename_prefix_of_fn(a_fa_fn), basename_prefix_of_fn(b_fa_fn)))
    b2a_m4_fn = op.join(work_dir, '%s.vs.%s.m4' % (basename_prefix_of_fn(b_fa_fn), basename_prefix_of_fn(a_fa_fn)))
    execute(blasr_cmd(query_fn=a_fa_fn, target_fn=b_fa_fn, out_fn=a2b_m4_fn))
    execute(blasr_cmd(query_fn=b_fa_fn, target_fn=a_fa_fn, out_fn=b2a_m4_fn))
    start, end = get_span_region_of_seq_from_m4(a2b_m4_fn=a2b_m4_fn, b2a_m4_fn=b2a_m4_fn, a_seq_name=a_seq_name)
    return (a_seq_name, start, end)


def get_span_region_of_seq_from_m4(a2b_m4_fn, b2a_m4_fn, a_seq_name):
    """
    a2b_m4_fn --- m4 output mapping FASTA file A to FASTA file B
    b2a_m4_fn --- m4 output mapping FASTA file B to FASTA file A
    a_seq_name --- name of a sequence (or chromosome) in FASTA file A
    return (a_seq_name, start, end) which spans all alignments in sequence of name `a_seq_name` in FASTA file A."""
    a2b_dict = _get_query_span_regions_from_m4(a2b_m4_fn)  # dict {query_name: NSE(query_name, start, end)}
    b2a_dict = _get_target_span_regions_from_m4(b2a_m4_fn)  # dict {target_name: NSE(target_name, start, end)}
    assert a_seq_name in a2b_dict.keys()
    assert a_seq_name in b2a_dict.keys()
    a2b_start, a2b_end = a2b_dict[a_seq_name].start, a2b_dict[a_seq_name].end
    b2a_start, b2a_end = b2a_dict[a_seq_name].start, b2a_dict[a_seq_name].end
    return apply_operator(a2b_start, b2a_start, min), apply_operator(a2b_end, b2a_end, max)


class NSE(object):
    def __init__(self, name=None, start=None, end=None):
        self.name = name
        self.start = start
        self.end = end


def _get_query_span_regions_from_m4(m4_fn):
    d = defaultdict(lambda: NSE())
    for r in BLASRM4Reader(m4_fn):
        if '/0_' not in r.qID:
            raise ValueError("M4 file %s must be generated by blasr align two FASTA files!" % (m4_fn))
        name = '/'.join(r.qID.split('/')[0:-1])  # BLASR adds '/0_len' to query fasta reads
        start, end = d[name].start, d[name].end
        start = r.qStart if start is None else min(start, r.qStart)
        end = r.qEnd if end is None else max(end, r.qEnd)
    return d


def _get_target_span_regions_from_m4(m4_fn):
    d = defaultdict(lambda: NSE())
    for r in BLASRM4Reader(m4_fn):
        name = r.sID
        start, end = d[name].start, d[name].end
        if r.sStrand == '+':
            start = r.sStart if start is None else min(start, r.sStart)
            end = r.sEnd if end is None else max(end, r.sEnd)
        elif r.sStrand == '-':
            start = (r.sLength - r.sEnd) if start is None else min(start,
                                                                   (r.sLength - r.sEnd))  # exclusive end to inclusive start
            end = (r.sLength - r.sStart) if end is None else max(end,
                                                                 (r.sLength - r.sStart))  # inclusive start to exclusive end
        d[name] = NSE(name, start, end)
    return d


def run(args):
    """Build consesus sequences from subreads.bam input using pbdagcon."""
    input_subreads_bam, ref_fa, output_prefix, consensus_id = args.input_subreads_bam, args.ref_fa, args.output_prefix, args.consensus_id,
    nproc, maxScore, use_first_seq_if_fail = args.nproc, args.maxScore, args.use_first_seq_if_fail

    sr_fasta = get_fasta_fn_from_subreads_bam_fn(input_subreads_bam)  # convert subreads.bam to subread.fasta
    pbdagcon_wrapper(fasta_filename=sr_fasta, output_prefix=output_prefix + '.all', consensus_name=consensus_id,
                     nproc=nproc, maxScore=maxScore, use_first_seq_if_fail=use_first_seq_if_fail)

    in_fa_fn, out_fa_fn = output_prefix + '.all.fasta', output_prefix + '.fasta'
    if ref_fa:  # use reference fasta to bound substr
        in_fa_obj = Fastafile(in_fa_fn)
        assert_fasta_has_one_seq(in_fa_obj)
        make_substr_fasta_of_seq_in_a_match_b(a_fa_obj=in_fa_obj, b_fa_obj=Fastafile(
            ref_fa), a_seq_name=in_fa_obj.references[0], out_fa_fn=out_fa_fn)
    else:
        execute('ln %s %s' % (in_fa_fn, out_fa_fn))

    try:  # OK to fail indexing fasta, this may happen if fasta is empty
        make_fai(out_fa_fn)
    except Exception:
        pass


def make_substr_fasta_of_seq_in_a_match_b(a_fa_obj, b_fa_obj, a_seq_name, out_fa_fn):
    """
    Get one substring of sequence a_seq_name in FASTA file a_fa_obj that match FASTA file b_fa_obj and write to out_fa_fn.
    if no such substring is found, write an empty file to out_fa_fn.
    """
    chrom, start, end = get_region_of_seq_in_a_match_b(
        a_fa_obj=a_fa_obj, b_fa_obj=b_fa_obj, a_seq_name=a_seq_name, work_dir=op.dirname(out_fa_fn))
    if start is not None and end is not None:
        substr_fasta(fileobj=a_fa_obj, chrom=chrom, start=start, end=end, out_fa_fn=out_fa_fn)
    else:  # otherwise, output is empty
        open(out_fa_fn, 'w').write('')


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))


if __name__ == "__main__":
    main()
