#!/usr/bin/env python
import os.path as op
import logging

from pbsv.independent.utils import execute, realpath
from pbsv.io.VcfIO import BedReader, BedWriter

from .independent import Constants as C
from .io import SVPolishFiles
from .utils import bed2prefix, write_to_bash_file, substr_fasta, get_ref_extension_for_sv

log = logging.getLogger()

collect_desc = 'Collect polished structural variants in directory.'


def run_collect(in_bed_fn, work_dir,
                polished_bed_fn, rejected_bed_fn, insufficient_bed_fn,
                min_qv, ref_ext_len):
    """
    Collect polished structural variants in work_dir.
    Write structural variants which do not have enough coverage, insufficient_coverage.bed.
    Write rejected structural variants to rejected.bed.
    Write polished structural variants to polished.bed.
    """
    reader = BedReader(in_bed_fn)
    p_writer = BedWriter(polished_bed_fn, samples=reader.samples)
    r_writer = BedWriter(rejected_bed_fn, samples=reader.samples)
    i_writer = BedWriter(insufficient_bed_fn, samples=reader.samples)

    for bed_record in reader:
        sv_prefix = bed2prefix(bed_record)
        log.info("Collecting polished SV for %s" % sv_prefix)
        data_dir = realpath(op.join(work_dir, sv_prefix))
        polished_bed_fn = SVPolishFiles(root_dir=data_dir, min_qv=min_qv, ref_ext_len=ref_ext_len).polish_ngmlr_bed
        if not op.exists(data_dir):
            log.info("No Polished structural variant detected, use the original one: %s " %
                     ' '.join(str(bed_record).split()[0:5]))
            i_writer.writeRecord(bed_record.to_str(reader.samples))
        else:
            polished_bed_records = [r for r in BedReader(polished_bed_fn)]
            if len(polished_bed_records) == 0:
                log.info("No Polished structural variant detected, use the original one: %s " %
                         ' '.join(str(bed_record).split()[0:5]))
                writer.writeRecord(bed_record.to_str(reader.samples))
            else:
                for r in polished_bed_records:
                    writer.writeRecord(r.to_str(reader.samples))

    reader.close()
    p_writer.close()
    r_writer.close()
    i_writer.close()
