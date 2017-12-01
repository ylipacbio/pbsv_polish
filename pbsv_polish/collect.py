#!/usr/bin/env python
from argparse import ArgumentParser
import os.path as op
import sys
import logging

from pbsv.independent.utils import execute, realpath
from pbsv.io.VcfIO import BedReader, BedWriter

from .independent import Constants as C
from .io import SVPolishFiles
from .utils import bed2prefix, write_to_bash_file, substr_fasta, get_ref_extension_for_sv

log = logging.getLogger()

collect_desc = 'Collect polished structural variants in directory.'


def run(args):
    run_collect(args.in_bed_fn,  args.out_dir, args.collected_bed_fn, args.min_qv, args.ref_ext_len)


def run_collect(in_bed_fn, out_dir, collected_bed_fn, min_qv, ref_ext_len):
    reader = BedReader(in_bed_fn)
    writer = BedWriter(collected_bed_fn, samples=reader.samples)

    for bed_record in reader:
        sv_prefix = bed2prefix(bed_record)
        log.info("Collecting polished SV for %s" % sv_prefix)
        data_dir = realpath(op.join(out_dir, sv_prefix))
        svp_files_obj = SVPolishFiles(root_dir=data_dir, min_qv=min_qv, ref_ext_len=ref_ext_len)
        polished_bed_fn = svp_files_obj.polish_ngmlr_bed

        if not op.exists(polished_bed_fn):
            log.info("No Polished structural variant detected, use the original one: %s " %
                  ' '.join(str(bed_record).split()[0:5]))
            writer.writeRecord(bed_record.to_str(reader.samples))
        else:
            polished_bed_records = [r for r in BedReader(polished_bed_fn)]
            if len(polished_bed_records) == 0:
                log.info("No Polished structural variant detected, use the original one: %s " %
                      ' '.join(str(bed_record).split()[0:5]))
                writer.writeRecord(bed_record.to_str(reader.samples))
            else:
                for r in polished_bed_records:
                    writer.writeRecord(r.to_str(reader.samples))

    writer.close()
