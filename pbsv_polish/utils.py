#!/usr/bin/env python

import os.path as op
import sys
import os
from collections import defaultdict
from pbcore.io import DataSet, FastaWriter

from pbsv.libs import Fastafile
from pbsv.io.linefile import X2PysamReader, iter_within_ref_regions
from pbsv.io.VcfIO import BedReader
from pbsv.independent.common import RefRegion
from pbsv.independent.utils import execute, realpath, write_to_bash_file, execute_as_bash, mv_cmd
from pbsv.cli import _mkdir


def f(strs):
    return '.'.join(strs)

#in_dir = 'in_sl_1813_yeast_10fold'

class SVInputFiles(object):
    def __init__(self, root_dir):
        def g(fn):
            return op.join(self.root_dir, fn)
        self.root_dir = realpath(root_dir)
        self.aln_bam = g('alignments.bam')
        self.subreads_xml = g("subreads.xml")
        self.genome_fa = g("genome.fa")
        self.sv_bed = g('structural_variants.bed')

class SVPolishFiles(object):
    def __init__(self, root_dir, min_qv=20):
        """root_dir is a directory containing all files for sv polishing."""
        def g(fn):
            return op.join(self.root_dir, fn)

        self.root_dir = realpath(root_dir)

        self.make_consensus_sh = g('make_consensus.sh')
        self.all_scripts_sh = g('all_scripts.sh')

        self._subreads_prefix = 'sr'
        self.subreads_bam = g(f([self._subreads_prefix, 'subreads.bam']))
        self.subreads_fa = g(f([self._subreads_prefix, 'subreads.fasta']))
        self.subreads_blasr_bam = g(f[self._subreads_prefix, 'blasr', 'bam'])

        self._sv_pbdagcon_prefix = 'sv_pbdagcon'
        self.sv_pbdagcon_fa = g(f([self._sv_pbdagcon_prefix, 'fasta']))

        self._sv_ref_prefix = 'sv_ref_w_ext'
        self.sv_ref_fa = g(f[self._sv_ref_prefx, 'fasta'])

        self._polish_prefix = 'polish.qv%s' % min_qv
        self.polish_fa = g(f([self._polish_prefix, 'fasta']))
        self.polish_fq = g(f([self._polish_prefix, 'fastq']))
        self.polish_blasr_bed = g(f([self._polish_prefix, 'blasr', 'bed']))
        self.polish_blasr_bam = g(f([self._polish_prefix, 'blasr', 'bam']))
        self.polish_blasr_sh = g(f([self._polish_prefix, 'blasr', 'sh']))

        self.polish_ngmlr_bed = g(f([self._polish_prefix, 'ngmlr', 'bed']))
        self.polish_ngmlr_sh = g(f([self._polish_prefix, 'ngmlr', 'sh']))
        self.polish_ngmlr_bam = g(f([self._polish_prefix, 'ngmlr', 'bam']))

        self._polish_hqlq_prefix = f(self._polish_prefix, 'hqlq')
        self.polish_hqlq_fa = g(f([self._polish_hqlq_prefix, 'fasta']))
        self.polish_hqlq_fq = g(f([self._polish_hqlq_prefix, 'fastq']))

def get_aln_reader(aln_fn, bed_fn):
    # reader = get_aln_reader(aln_fn=aln_fn, bed_fn=bed_fn)
    ref_regions = get_ref_regions_from_bed_file(bed_fn)
    reader = X2PysamReader(aln_fn, ref_regions)
    return reader

def yield_alns_from_bed_file(alnfile_obj, bedreader_obj):
    for bed_record in bedreader_obj:
        ref_region = to_ref_region(bed_record)
        yield (bed_record, get_alns_within_ref_region(alnfile_obj, ref_region))

def to_ref_region(bed_record):
    return RefRegion(bed_record.chrom, bed_record.start, bed_record.end+1)

def get_ref_regions_from_bed_file(bed_fn):
    return [to_ref_region(bed_r) for bed_r in BedReader(bed_fn)]

def yield_ref_region_from_bed_file(bed_fn):
    for bed_r in BedReader(bed_fn):
        yield to_ref_region(bed_r)

def get_alns_within_ref_region(alnfile_obj, ref_region):
    """Return a list of alignments within a ref region."""
    return [aln for aln in iter_within_ref_regions(alnfile_obj, [ref_region])]

def zmw_from_subread(subread):
    """Given a subread 'movie/zmw/start_end', return 'movie/zmw'"""
    return '/'.join(subread.split('/')[0:2])

def get_query_zmws_from_alns(alns):
    """Given a list of alignments, return a list of non-redundant query zmws"""
    return list(set([zmw_from_subread(sr) for sr in get_query_subreads_from_alns(alns)]))

def get_query_subreads_from_alns(alns):
    """Given a list of alignments, return a list of non-redundant query subreads"""
    return list(set([aln.query_name for aln in alns]))

from svkits.utils import get_movie2zmws_from_zmws, make_subreads_bam, make_subreads_bam_using_pbcore

def make_subreads_bam_of_zmws(movie2bams, zmws, out_prefix, dry_run=False):
    movie2zmws = get_movie2zmws_from_zmws(zmws)
    return make_subreads_bam(movie2zmws, movie2bams, out_prefix, dry_run=dry_run)
    #TODO
    #return make_subreads_bam_using_pbcore(movie2zmws, movie2bams, out_prefix, dry_run=dry_run)

def get_subreads_bam_files_from_xml(in_subreads_xml):
    return [fn for fn in DataSet(in_subreads_xml).toExternalFiles() if fn.endswith('subreads.bam')]

def get_movies2bams_from_subreads_xml(in_subreads_xml):
    subreads_fns = get_subreads_bam_files_from_xml(in_subreads_xml)
    return get_movies2bams_from_subreads_bam_files(subreads_fns)

def get_movies2bams_from_subreads_bam_files(subreads_bam_fns):
    """Input: a list of subreads.bam files.
       Return {movie: bam_fn}
    """
    movie2bams = defaultdict(lambda:set())
    for fn in subreads_bam_fns:
        movie = fn.split('/')[-1].split('.')[0]
        print "movie=%s" % movie
        if movie in movie2bams and not fn == movie2bams[movie]:
            raise ValueError("Movie %s mapping to multiple bam files %r and %r" % (movie, movie2bams[movie], fn))
        movie2bams[movie].add(fn)
    return movie2bams

def sort_index_bam_inline(in_fn, out_fn, nproc=4, tmp_dir=None):
    """Sort and make index of sam or bam files inline"""
    sort_arg = '-T {tmp_dir}/srt' if tmp_dir is not None else ''
    tmp_fn = out_fn + '.sorted'
    c0 = "samtools sort --threads {nproc} {sort_arg} {in_fn} -o {tmp_fn}".format(
          nproc=nproc, sort_arg=sort_arg, in_fn=in_fn, tmp_fn=tmp_fn)
    c1 = mv_cmd(tmp_fn, out_fn)
    c2 = "samtools index {out_fn}".format(out_fn=out_fn)
    return " ; ".join([c0, c1, c2])
