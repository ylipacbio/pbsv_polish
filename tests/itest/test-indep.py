from nose.tools import (assert_equal, assert_false, assert_true, assert_raises)
from pbsv.independent.common import (_SvEnumCls, SvFmt, SvType, SvAnnot, SvAnnotations, _2gt, _f, rev_dict,
        ReadSvWindow, ReadSvWindowCluster, RefRegion, to_ref_regions, get_chrom_start_end_from_string)
from pbsv.independent.annot import Repeat, annot_seq, REPEATS
from StringIO import StringIO
import pbsv.independent.basic as basic
import pbsv.independent.iobase as iobase
import pbsv.independent.utils as utils
import pbsv.independent.FastaReader as FastaReader
import pytest
import unittest
import collections
import gzip
import json
import os
import os.path as op

THIS_DIR = op.dirname(os.path.abspath(os.path.realpath(__file__)))
DATA_DIR = op.join(THIS_DIR, '../data')
OUT_DIR = op.join(THIS_DIR, '../out')
PBSV_POLISH_CFG = op.join(DATA_DIR, 'pbsv.polish.cfg')

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)
log = logging.getLogger()

if not os.path.isdir(OUT_DIR):
    os.makedirs(OUT_DIR)

def test_iobase_getFileHandle():
    fn = op.relpath(op.join(OUT_DIR, 'test_iobase_getFileHandle.gz'))
    assert isinstance(iobase.getFileHandle(fn, mode='w'), gzip.GzipFile)
    with pytest.raises(Exception) as excinfo:
        iobase.getFileHandle(object())
    assert "Invalid type <type 'object'> to getFileHandle" == str(excinfo.value)

from pbsv_polish.utils import sort_index_bam_inline_cmd, blasr_cmd, _fn2fmtarg,\
        variant_caller_cmd, trim_lq_cmd, bed2prefix
from pbsv_polish.utils import *

def test_sort_index_bam_inline_cmd():
    s = sort_index_bam_inline_cmd('1.bam')
    assert s == 'samtools sort --threads 4  1.bam -o 1.bam.sorted ; mv 1.bam.sorted 1.bam ; samtools index 1.bam'

def test_fn2fmtarg():
    assert _fn2fmtarg('out.m4') == '-m 4'
    assert _fn2fmtarg('out.m0') == '-m 0'
    assert _fn2fmtarg('out.bam') == '--bam'

def test_blasr_cmd():
    s = blasr_cmd('q.fa', 't.fa', 'out.bam')
    assert s == 'blasr q.fa t.fa --bam --out out.bam --nproc 8 --maxMatch 15 --bestn 10 --hitPolicy randombest'

def test_variant_caller_cmd():
    s = variant_caller_cmd(align_bam='a.bam', ref_fa='ref.fa', nproc=8, out_fa='o.fa', out_fq='o.fq')
    assert s ==  'variantCaller --algorithm=best a.bam --verbose -j8 --referenceFilename=ref.fa -o o.fa -o o.fq'

def test_trim_lq_cmd():
    s = trim_lq_cmd('i.fq', 'o.fq', 'o.fa', 30)
    assert s == 'trim_lq i.fq o.fq --min_qv 30 ; fq2fa o.fq o.fa'

def test_bed2prefix():
    from pbsv.io.VcfIO import BedRecord
    o = BedRecord(chrom='chr1', start=0, end=100, sv_type='Deletion',
                  sv_len=-100, seq=None, fmt='0/1:3:6', annotations=['ALU'])
    s = bed2prefix(o)
    assert s == 'chr1_0_100_Deletion_-100'

def test_sv_transform_coordinate_cmd():
    s = sv_transform_coordinate_cmd('i.bed', 'o.bed')
    assert s == 'sv_transform_coordinate i.bed o.bed'

def test_pbsv_align_cmd():
    s = pbsv_align_cmd(reads_fn='q.bam', ref_fa_fn='ref.fa', cfg_fn='my.cfg', o_bam_fn='o.bam')
    assert s == 'pbsv align ref.fa q.bam o.bam --cfg_fn=my.cfg'

def test_pbsv_run_cmds():
    s = pbsv_run_and_transform_cmds('read.bam', 'ref.fa', PBSV_POLISH_CFG, 'o.bam', 'o.bed', 'ngmlr')
    assert s[0].startswith('pbsv align')
    assert s[1].startswith('pbsv call')
    assert s[2].startswith('sv_transform_coordinate')
    assert s[3].startswith('samtools sort')

    s = pbsv_run_and_transform_cmds('read.bam', 'ref.fa', PBSV_POLISH_CFG, 'o.bam', 'o.bed', 'blasr')
    assert s[0].startswith('blasr')
    assert s[1].startswith('pbsv call')
    assert s[2].startswith('sv_transform_coordinate')
    assert s[3].startswith('samtools sort')

def test_SVPolishFiles():
    s = SVPolishFiles('dir')
    assert s.subreads_bam == 'dir/sr.subreads.bam'
    assert s.subreads_consensus_sh == 'dir/subreads_consensus.sh'
    assert s.dagcon_fa == 'dir/sv_pbdagcon.fasta'
