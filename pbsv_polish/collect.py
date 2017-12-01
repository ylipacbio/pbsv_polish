#!/usr/bin/env python
import os.path as op
import logging

from pbsv.independent.utils import execute, realpath
from pbsv.io.VcfIO import BedReader, BedWriter

from .independent import Constants as C
from .io import SVPolishFiles
from .utils import bed2prefix, write_to_bash_file, substr_fasta, get_ref_extension_for_sv, ActionRecord

log = logging.getLogger()

collect_desc = 'Collect polished structural variants in directory.'


def run_collect(work_dir, collected_bed_fn, min_qv, ref_ext_len):
    """
    Collect polished structural variants in work_dir, while work_dir must contain
    * in.passed.bed --- BED file containing original input structural variants
        which have passed min coverage filter
    * in.skipped.bed --- BED file containing original input structural variants
        which were skipped because not having enough coverage
    Then `run_collect` will write structural variants into different BED files by categories:
    * in.rejected.bed --- BED file containing original input structural variants
        which are rejected (e.g., not having consensus sequence) by polishing
    * in.polished.bed --- BED file containing original input structural variants
        which have been polished
    * out.polished.bed --- BED file containing polished structural variants
    """
    in_passed_bed_fn, in_skipped_bed_fn= op.join(work_dir, 'in.passed.bed'), op.join(work_dir, 'in.skipped.bed')
    in_rejected_bed_fn =  op.join(work_dir, 'in.rejected.bed')
    in_polished_bed_fn = op.join(work_dir, 'in.polished.bed')
    out_polished_bed_fn = op.join(work_dir, 'out.polished.bed')

    reader = BedReader(in_passed_bed_fn)
    in_polished_writer = BedWriter(in_polished_bed_fn, samples=reader.samples)
    in_rejected_writer = BedWriter(in_rejected_bed_fn, samples=reader.samples)
    out_polished_writer = BedWriter(out_polished_bed_fn, samples=[C.CONSENSUS_SAMPLE])
    collected_writer = BedWriter(collected_bed_fn, samples=[C.CONSENSUS_SAMPLE])

    prev_svobj = None

    for bed_record in reader:
        sv_prefix = bed2prefix(bed_record)
        svobj_dir = realpath(op.join(work_dir, sv_prefix))
        polished_bed_fn = SVPolishFiles(root_dir=svobj_dir, min_qv=min_qv, ref_ext_len=ref_ext_len).polish_ngmlr_bed
        polished_bed_records = [r for r in BedReader(polished_bed_fn)]
        if len(polished_bed_records) == 0:
            log.info(str(ActionRecord.from_obj(bed_record, action=ActionRecord.REJECTED, comment='NO_SV')))
            in_rejected_writer.writeRecord(bed_record)
        else:
            log.info(str(ActionRecord.from_obj(bed_record, action=ActionRecord.POLISHED, comment='POLISHED')))
            in_polished_writer.writeRecord(bed_record)
            for r in polished_bed_records:
                out_polished_writer.writeRecord(r)
                if prev_svobj and is_similar_svobj(r, prev_svobj):
                    log.info('Merge {} with {}'.format(bed2prefix(r), bed2prefix(prev_svobj)))
                else:
                    collected_writer.writeRecord(r)
                prev_svobj = r

    reader.close()
    in_polished_writer.close()
    in_rejected_writer.close()
    out_polished_writer.close()
    collected_writer.close()

    def _f(bed_fn, what):
        n = len([r for r in BedReader(bed_fn)])
        log.info("There are {n} {what} structural variants in file {f} .".format(n=n, what=what, f=bed_fn))

    _f(in_passed_bed_fn, 'PASSED')
    _f(in_skipped_bed_fn, 'SKIPPED')
    _f(in_rejected_bed_fn, 'REJECTED')
    _f(in_polished_bed_fn, 'INPUT POLISHED')
    _f(out_polished_bed_fn, 'OUTPUT POLISHED')
    _f(collected_bed_fn, 'COLLECTED POLISHED')
