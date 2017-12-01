#!/usr/bin/env python

import os.path as op
import numpy as np

from pbcore.io import FastaWriter, SubreadSet
from pbsv.libs import Fastafile, AlignmentFile
from pbsv.io.bamstream import iterator_of_alignments_in_ref_regions
from pbsv.io.VcfIO import BedRecord
from pbsv.independent.common import RefRegion, SvType, SvFmts, SvFmt
from pbsv.independent.utils import _is_fmt, cmds_to_bash, execute, realpath, autofmt, is_fasta
from .independent import Constants as C

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()


def get_movie_and_zmw_from_name(name):
    """Given a string of pacbio zmw name or read name, return movie and zmw"""
    try:
        fs = name.strip().split(' ')[0].split('/')
        movie, zmw = fs[0], fs[1]
        return movie, int(zmw)
    except ValueError:
        raise ValueError("Read %r is not a PacBio read." % name)


def write_fasta(out_fa_fn, records):
    """Write a list of fasta records [(name, seq), ...,  (name, seq)] to out_fa_fn"""
    with FastaWriter(out_fa_fn) as w:
        for r in records:
            w.writeRecord(r[0], r[1])


def substr_fasta(fileobj, chrom, start, end, out_fa_fn):
    """fetch a substring of reference fasta sequence and save to output fasta file out_fa_fn"""
    try:
        seq = fileobj.fetch(str(chrom), int(start), int(end))
    except Exception:
        raise ValueError("Could not get substring (%s, %s, %s) from %s" % (chrom, start, end, fileobj.filename))
    name = '%s__substr__%s_%s' % (chrom, start, end)
    write_fasta(out_fa_fn, [(name, seq)])


def yield_alns_from_bed_file(alnfile_obj, bedreader_obj):
    for bed_record in bedreader_obj:
        ref_region = to_ref_region(bed_record)
        yield (bed_record, get_alns_within_ref_region(alnfile_obj, ref_region))


def to_ref_region(bed_record):
    return RefRegion(bed_record.chrom, bed_record.start, bed_record.end + 1)


def get_alns_within_ref_region(alnfile_obj, ref_region):
    """Return a list of alignments within a ref region."""
    return [aln for aln in iterator_of_alignments_in_ref_regions(alnfile_obj, [ref_region])]


def zmw_from_subread(subread):
    """Given a subread 'movie/zmw/start_end', return 'movie/zmw'"""
    try:
        return '/'.join(subread.split('/')[0:2])
    except Exception:
        raise ValueError("Could not convert read %s to zmw" % subread)


def write_to_bash_file(cmds, bash_sh_fn, write_mode='w'):
    """Write commands to a bash script."""
    with open(bash_sh_fn, write_mode) as writer:
        writer.write(cmds_to_bash(cmds) + '\n')


def get_query_zmws_from_alns(alns):
    """Given a list of alignments, return a list of non-redundant query zmws"""
    return list(set([zmw_from_subread(sr) for sr in get_query_subreads_from_alns(alns)]))


def get_query_subreads_from_alns(alns):
    """Given a list of alignments, return a list of non-redundant query subreads"""
    return list(set([aln.query_name for aln in alns]))


def sort_zmws_by_moive(zmws):
    """ sort a list of zmws by moive names."""
    return sorted(zmws, key=lambda zmw: get_movie_and_zmw_from_name(zmw)[0])


def _get_subreads_ds_from_fn_or_obj(in_subreads_fn_or_obj):
    return SubreadSet(in_subreads_fn_or_obj) if isinstance(in_subreads_fn_or_obj, str) else in_subreads_fn_or_obj


def get_fastafile_obj_from_fn_or_obj(in_fa_fn_or_or_obj):
    return Fastafile(in_fa_fn_or_or_obj) if isinstance(in_fa_fn_or_or_obj, str) else in_fa_fn_or_or_obj


def get_bam_header_from_subreads_ds(ds):
    """Return a BamHeader obj from a SubreadSet.
    ds --- SubreadSet
    """
    from pbtranscript.io.PbiBamIO import BamHeader
    from pbcore.io import IndexedBamReader
    header = BamHeader(ignore_pg=True)
    for bam in ds.resourceReaders():
        if not isinstance(bam, IndexedBamReader):
            raise ValueError("%s must have pbi index generated", bam.filename)
        header.add(bam.peer.header)
        for rg in bam.peer.header["RG"]:  # readGroupTable:
            if rg['PL'] != "PACBIO":
                raise IOError("Input BAM file %s for %s must be PacBio BAM.", bam.filename, 'get_bam_header')
        for rg in bam.readGroupTable:
            assert rg.ReadType in ["SUBREAD"]
    return header.header  # a dict


def make_subreads_bam_of_zmws2(in_subreads_fn_or_obj, zmws, out_bam_fn, out_fa_fn=None):
    """Make subreads bam of zmws
    in_subreads_fn --- subreads.bam or suberadset.xml
    # movie2bams --- {movie: bam_file}, e.g., {'movie1': 'path-to-movie1', 'movie2': 'path-to-movie2'}
    zmws --- a list of zmws or reads, e.g. ['movie1/100', 'movie2/200'], or ['movie1/100/0_100']
    if out_fa_fn is not None, generate fasta file as well
    """
    subreads_ds = _get_subreads_ds_from_fn_or_obj(in_subreads_fn_or_obj)
    assert isinstance(subreads_ds, SubreadSet)
    zmws = sort_zmws_by_moive(zmws)
    header = get_bam_header_from_subreads_ds(subreads_ds)
    reads = [r for r in yield_subreads_of_zmws_in_ds(subreads_ds=subreads_ds, zmws_or_reads=zmws)]
    write_reads_bam(out_bam_fn=out_bam_fn, header=header, reads=reads)
    if out_fa_fn:
        write_reads_fasta(out_fa_fn=out_fa_fn, reads=reads)


def get_non_redundant_zmws_from_zmws_or_reads(zmws_or_reads):
    """Return a list of non-redundant zmws from a list of zmws or reads"""
    return sort_zmws_by_moive(list(set(['/'.join([str(x) for x in get_movie_and_zmw_from_name(zmw)]) for zmw in zmws_or_reads])))


def yield_subreads_of_zmws_in_ds(subreads_ds, zmws_or_reads):
    """
    Iterate over subreads of zmws in a subreads dataset
    subreads_ds -- a SubreadSet object
    zmws_or_reads --- a list of zmws or reads, e.g., ['movie1/100', 'movie2/200'] or ['movie1/100/0_10']
    """
    assert isinstance(subreads_ds, SubreadSet)
    zmws = get_non_redundant_zmws_from_zmws_or_reads(zmws_or_reads)
    for zmw in zmws:
        for subread in subreads_of_a_zmw_in_ds(subreads_ds, zmw):
            yield subread


def subreads_of_a_zmw_in_ds(subreads_ds, zmw):
    """
    Return a list of all subreads of a zmw in subreads_ds
    subreads_ds -- a SubreadSet object
    zmws --- a list of zmws, e.g., e.g. ['movie1/100', 'movie2/200']
    """
    movie, zmw_int = get_movie_and_zmw_from_name(zmw)
    rows = np.nonzero(np.logical_and(  # pylint: disable=no-member
        subreads_ds.index.qId == subreads_ds.movieIds[movie],
        subreads_ds.index.holeNumber == zmw_int))[0]
    return subreads_ds[rows]


def bed2prefix(bed_record):
    fields = [bed_record.chrom, bed_record.start, bed_record.end, bed_record.sv_type, bed_record.sv_len]
    return '_'.join([str(x) for x in fields])


def qsub_to_sge_or_run_local(script_fn, use_sge=False):
    if use_sge:
        cmd = 'qsub -q default -V -pe smp 1 -v -cwd -b y {f} -e {f}.qerr.log -o {f}.qout.log 1>{f}.stdout.log 2>{f}.stderr.log'.format(
            f=script_fn)
        log.info(cmd)
    else:
        cmd = 'chmod +x %s && bash %s' % (script_fn, script_fn)
        log.info(cmd)
    execute(cmd)


def write_reads_bam(out_bam_fn, header, reads):
    with AlignmentFile(out_bam_fn, "wb", header=header, check_sq=False) as writer:
        for subread in reads:
            writer.write(subread.peer)


def write_reads_fasta(out_fa_fn, reads):
    with FastaWriter(out_fa_fn) as writer:
        for read in reads:
            writer.writeRecord(read.readName, read.read(aligned=False))


def get_ref_extension_for_sv(bed_record, ref_seq_len=None):
    """Get reference extension for structural variant"""
    s = max(0, bed_record.start - max(C.REFERENCE_EXTENSION_SV_FACTOR * bed_record.sv_len, C.REFERENCE_EXTENSION))
    e = bed_record.end + max(C.REFERENCE_EXTENSION_SV_FACTOR * bed_record.sv_len, C.REFERENCE_EXTENSION)
    if ref_seq_len is not None:
        e = min(e, int(ref_seq_len))
    return s, e


def assert_fasta_has_one_seq(fasta_obj_or_fn):
    obj = get_fastafile_obj_from_fn_or_obj(fasta_obj_or_fn)
    if len(obj.references) != 1:
        raise ValueError("FASTA file %s must contain exactly one sequence, while it contains: %s" %
                         (obj.filename, '\n'.join(obj.references)))


def apply_operator(a, b, op):
    if (a is None and b is None):
        return None
    elif a is None:
        return b
    elif b is None:
        return a
    else:
        return op(a, b)


def basename_prefix_of_fn(fn):
    """Return base name prefix of a file path, e.g.,
    ..doctest:
        >>> basename_prefix_of_fn('/home/my.b.txt')
        'my.b'
    """
    basename = op.basename(fn)
    return basename[0:basename.rfind('.')] if '.' in basename else basename


def get_sv_from_non_pbsv_bed(in_bed):
    """Read in_bed which has a different format from pbsv BED as if it is a pbsv BED.
    Return a list of BedRecord objects.
    Example of in_bed:
    chr1\t1000\t1000\tInsertion\t100\t.\t.\t.
    Each line must contain at least five columns: chromsome, start, end, svtype, svlen
    """
    ret = []
    for idx, line in enumerate(open(in_bed, 'r')):
        fs = line.strip().split('\t')
        try:
            chr, start, end, svtype, svlen = fs[0], int(fs[1]), int(fs[2]), fs[3], int(fs[4])
            if SvType(svtype).is_Deletion:
                svlen = -abs(svlen)
        except:
            raise AssertionError(
                'True set BED line {idx} must contain at least five columns, (chr, start, end, type, len),          line={line}'.format(idx=idx, line=line))
        fmts = SvFmts([SvFmt.fromString('0/0:1:1')], ['SAMPLE'])
        ret.append(BedRecord(chrom=chr, start=start, end=end, sv_id='.', sv_type=svtype,
                             sv_len=svlen, alt=None, annotations=None,   fmts=fmts))
    return ret


class ActionRecord(object):
    """
    Simple action record class to record actions taken for a structural variant
    """
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

