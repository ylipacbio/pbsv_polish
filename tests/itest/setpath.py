#!/usr/bin/env python
import os.path as op

THIS_DIR = op.dirname(op.abspath(__file__))
ROOT_DIR = op.dirname(THIS_DIR)
OUT_DIR = op.join(ROOT_DIR, 'out')
DATA_DIR = op.join(ROOT_DIR, 'data')

SIV_DIR = '/pbi/dept/secondary/siv/testdata/pbsvp-unittest/'
SIV_DATA_DIR = op.join(SIV_DIR, 'data')
SIV_EXPECTED_DIR = op.join(SIV_DIR, 'expected_out') # new test stdout
