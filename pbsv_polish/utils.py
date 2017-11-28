#!/usr/bin/env python

import os.path as op
from collections import defaultdict
import numpy as np

from pbcore.io import DataSet, FastaWriter, SubreadSet

from pbsv.libs import Fastafile, AlignmentFile
from pbsv.io.bamstream import BamStream, iterator_of_alignments_in_ref_regions
from pbsv.io.VcfIO import BedReader, BedRecord
from pbsv.independent.common import RefRegion, SvType, SvFmts, SvFmt
from pbsv.independent.utils import _is_fmt, cmds_to_bash, execute, realpath, mv_cmd, autofmt, is_fasta
from pbsv.run import svcall_cmd, ngmlrmap_cmd
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


def is_fastq(fn):
    """Return true if a file extension is fq or fastq"""
    return _is_fmt(fn, ["fq", "fastq"])


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


def get_aln_reader(aln_fn, bed_fn):
    ref_regions = get_ref_regions_from_bed_file(bed_fn)
    reader = BamStream(fn=aln_fn, ref_regions=ref_regions, require_sorted=True)
    return reader


def yield_alns_from_bed_file(alnfile_obj, bedreader_obj):
    for bed_record in bedreader_obj:
        ref_region = to_ref_region(bed_record)
        yield (bed_record, get_alns_within_ref_region(alnfile_obj, ref_region))


def to_ref_region(bed_record):
    return RefRegion(bed_record.chrom, bed_record.start, bed_record.end + 1)


def get_ref_regions_from_bed_file(bed_fn):
    return [to_ref_region(bed_r) for bed_r in BedReader(bed_fn)]


def yield_ref_region_from_bed_file(bed_fn):
    for bed_r in BedReader(bed_fn):
        yield to_ref_region(bed_r)


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
    rows = np.nonzero(np.logical_and(subreads_ds.index.qId == subreads_ds.movieIds[movie], subreads_ds.index.holeNumber == zmw_int))[
        0]  # pylint: disable=no-member
    return subreads_ds[rows]


def get_subreads_bam_files_from_xml(in_subreads_xml):
    return [fn for fn in DataSet(in_subreads_xml).toExternalFiles() if fn.endswith('subreads.bam')]


def get_movies2bams_from_subreads_xml(in_subreads_xml):
    subreads_fns = get_subreads_bam_files_from_xml(in_subreads_xml)
    return get_movies2bams_from_subreads_bam_files(subreads_fns)


def get_movies2bams_from_subreads_bam_files(subreads_bam_fns):
    """Input: a list of subreads.bam files.
       Return {movie: bam_fn}
    """
    movie2bams = defaultdict(lambda: set())
    for fn in subreads_bam_fns:
        movie = fn.split('/')[-1].split('.')[0]
        if movie in movie2bams and not fn == movie2bams[movie]:
            raise ValueError("Movie %s mapping to multiple bam files %r and %r" % (movie, movie2bams[movie], fn))
        movie2bams[movie].add(fn)
    return movie2bams


def bed2prefix(bed_record):
    fields = [bed_record.chrom, bed_record.start, bed_record.end, bed_record.sv_type, bed_record.sv_len]
    return '_'.join([str(x) for x in fields])


def sort_index_bam_inline_cmd(in_fn, nproc=4, tmp_dir=None):
    """Sort and make index of sam or bam files inline"""
    sort_arg = '-T {tmp_dir}/srt' if tmp_dir is not None else ''
    tmp_fn = in_fn + '.sorted'
    c0 = "samtools sort --threads {nproc} {sort_arg} {in_fn} -o {tmp_fn}".format(
        nproc=nproc, sort_arg=sort_arg, in_fn=in_fn, tmp_fn=tmp_fn)
    c1 = mv_cmd(tmp_fn, in_fn)
    c2 = "samtools index {in_fn}".format(in_fn=in_fn)
    return " ; ".join([c0, c1, c2])


def pbindex_cmd(fn):
    return 'pbindex {fn}'.format(fn=fn)


def _fn2fmtarg(fn):
    fnext2fmtarg = {"m0": "-m 0", "m4": "-m 4", "bam": "--bam"}
    return fnext2fmtarg[autofmt(fn, fnext2fmtarg.keys())[1]]


def blasr_cmd(query_fn, target_fn, out_fn, nproc=8):
    return "blasr {q} {t} {fmt} --out {out_fn} --nproc {nproc} --maxMatch 15 --bestn 10 --hitPolicy randombest".\
        format(q=query_fn, t=target_fn, out_fn=out_fn, fmt=_fn2fmtarg(out_fn), nproc=nproc)


def sv_pbdagcon_cmd(subreads_bam, output_prefix, output_seq_id, ref_fa):
    return 'sv_pbdagcon %s %s %s --ref_fa %s' % (subreads_bam, output_prefix, output_seq_id, ref_fa)


def variant_caller_cmd(align_bam, ref_fa, out_fa, out_fq, nproc):
    return "variantCaller --algorithm=best {align_bam} --verbose -j{nproc} --referenceFilename={ref_fa} -o {out_fa} -o {out_fq}".\
        format(align_bam=align_bam, ref_fa=ref_fa, out_fa=out_fa, out_fq=out_fq, nproc=nproc)


def trim_lq_cmd(in_fq, out_fq, out_fa, min_qv):
    assert is_fastq(in_fq) and is_fastq(out_fq) and is_fasta(out_fa)
    # simply remove lower case sequences on both ends
    c0 = 'trim_lq {in_fq} {out_fq} --min_qv {min_qv}'.format(in_fq=in_fq, out_fq=out_fq, min_qv=min_qv)
    c1 = 'fq2fa {fq} {fa}'.format(fq=out_fq, fa=out_fa)
    return ' ; '.join([c0, c1])


def sv_transform_coordinate_cmd(in_sv_fn, o_sv_fn):
    return 'sv_transform_coordinate %s %s' % (in_sv_fn, o_sv_fn)


def pbsv_align_cmd(reads_fn, ref_fa_fn, cfg_fn, o_bam_fn):
    return 'pbsv align {r} {q} {o} --cfg_fn={c}'.format(r=ref_fa_fn, q=reads_fn, c=cfg_fn, o=o_bam_fn)


def pbsv_run_and_transform_cmds(reads_fn, ref_fa_fn, cfg_fn, o_bam_fn, o_bed_fn, algorithm='ngmlr'):
    """Using blasr to align reads to ref_fa_fn, then run `pbsv call` to call structural variants."""
    #c0 = _align_cmd(reads_fn=reads_fn, ref_fa_fn=ref_fa_fn, cfg_fn=cfg_fn, algorithm=algorithm)
    if algorithm == 'ngmlr':
        c0 = pbsv_align_cmd(reads_fn=reads_fn, ref_fa_fn=ref_fa_fn, cfg_fn=cfg_fn, o_bam_fn=o_bam_fn)
    elif algorithm == 'blasr':
        c0 = blasr_cmd(query_fn=reads_fn, target_fn=ref_fa_fn, out_fn=o_bam_fn)
    else:
        raise ValueError('Could not use algorithm %s to align reads to reference.' % algorithm)

    tmp_bed = o_bed_fn + '.use_substr_as_chrom.bed'
    c1 = svcall_cmd(ref_fn=ref_fa_fn, in_bam=o_bam_fn, out_bed=tmp_bed, cfg=cfg_fn)
    c2 = rm_ngmlr_indices_cmd(ref_fn=ref_fa_fn)
    c3 = sv_transform_coordinate_cmd(tmp_bed, o_bed_fn)
    c4 = sort_index_bam_inline_cmd(o_bam_fn)
    return [c0, c1, c2, c3, c4]


def rm_ngmlr_indices_cmd(ref_fn):
    from pbsv.ngmlrmap import _ngm_1, _ngm_2
    return 'rm -f %s %s' % (_ngm_1(ref_fn), _ngm_2(ref_fn))


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
        "my.b"
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
