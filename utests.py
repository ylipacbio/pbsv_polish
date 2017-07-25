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

from pbsv_polish.utils import sort_index_bam_inline
def test_sort_index_bam_inline():
    s = sort_index_bam_inline('1.bam', '2.bam')
    assert s == 'samtools sort --threads 4  1.bam -o 2.bam.sorted ; mv 2.bam.sorted 2.bam ; samtools index 2.bam'
