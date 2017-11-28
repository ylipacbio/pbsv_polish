#!/usr/bin/env python
import os.path as op


class Constants(object):
    PBSV_POLISH_CFG = op.join(op.dirname(op.dirname(__file__)), 'data', 'pbsv.polish.cfg')
    REFERENCE_EXTENSION = 1000
    REFERENCE_EXTENSION_SV_FACTOR = 2
    MIN_POLISH_QV = 20
    MIN_POLISH_COVERAGE = 5
    BLASR_NPROC = 4
    VARIANT_CALLER_NPROC = 4
