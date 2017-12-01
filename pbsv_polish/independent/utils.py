#!/usr/bin/env python

import os
import os.path as op
from pbsv.independent.utils import _is_fmt, cmds_to_bash, execute, realpath, autofmt, is_fasta


def is_fastq(fn):
    """Return true if a file extension is fq or fastq"""
    return _is_fmt(fn, ["fq", "fastq"])
