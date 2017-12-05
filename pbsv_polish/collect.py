#!/usr/bin/env python
import os.path as op
import logging

from pbsv.independent.utils import execute, realpath
from pbsv.io.VcfIO import BedReader, BedRecord, BedWriter

from .independent import Constants as C
from .io import SVPolishFiles
from .utils import bed2prefix, write_to_bash_file, substr_fasta, get_ref_extension_for_sv, ActionRecord

log = logging.getLogger()

collect_desc = 'Collect polished structural variants in directory.'
MAX_START_DIFF = 10 # maximum start position distance to merge structural variants
MAX_SVLEN_DIFF_PCT = 0.05 # maximum sv length percentage difference to merge structural variants

def cmp_bedrecord_f(aobj, bobj, max_start_diff, max_svlen_diff_pct):
    assert isinstance(aobj, BedRecord) and isinstance(bobj, BedRecord)
    if aobj.chrom == bobj.chrom:
        if abs(aobj.start - bobj.start) <= max_start_diff and \
            abs(aobj.sv_len - bobj.sv_len) * 1.0 / max(aobj.sv_len, bobj.sv_len) <= max_svlen_diff_pct:
            return 0
        return aobj.start - bobj.start
    else:
        return 1 if aobj.chrom > bobj.chrom else -1


def remove_redundant_objs(objs, cmp_f):
    """Remove redundant objects from objs and return non-redundant objects.
    objs --- a list of sorted objects.
    cmp_f --- cmp_f(aobj, bobj) return 0 if aobj 'equals' bobj and -1 if aobj is less than bobj, and 1 otherwise.
    """
    ret = []
    for obj in objs:
        if ret:
            if cmp_f(obj, ret[-1]) == 0:
                log.info('Redundant objects {} and {}'.format(repr(obj), repr(ret[-1])))
                continue
        ret.append(obj)
    return ret


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
    * out.polished.bed --- BED file containing polished structural variants, may
      contain redundant structural variants
    The output file collected_bed_fn is identical to out.polished.bed except that
    redundant structural variants have been removed.
    """
    in_passed_bed_fn, in_skipped_bed_fn = op.join(work_dir, 'in.passed.bed'), op.join(work_dir, 'in.skipped.bed')
    in_rejected_bed_fn = op.join(work_dir, 'in.rejected.bed')
    in_polished_bed_fn = op.join(work_dir, 'in.polished.bed')
    out_polished_bed_fn = op.join(work_dir, 'out.polished.bed')

    reader = BedReader(in_passed_bed_fn)
    in_polished_writer = BedWriter(in_polished_bed_fn, samples=reader.samples)
    in_rejected_writer = BedWriter(in_rejected_bed_fn, samples=reader.samples)
    out_polished_writer = BedWriter(out_polished_bed_fn, samples=[C.CONSENSUS_SAMPLE])

    collected_bed_records = []
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
            collected_bed_records.extend(polished_bed_records)

    def bed_record_to_key(bed_record):
        return (bed_record.chrom, bed_record.start)
    def cmp_f(aobj, bobj):
        return cmp_bedrecord_f(aobj, bobj, MAX_START_DIFF, MAX_SVLEN_DIFF_PCT)
    collected_bed_records = sorted(collected_bed_records, key=bed_record_to_key)
    collected_bed_records = remove_redundant_objs(collected_bed_records, cmp_f)
    with BedWriter(collected_bed_fn, samples=[C.CONSENSUS_SAMPLE]) as writer:
        for bed_record in collected_bed_records:
            writer.writeRecord(bed_record)

    reader.close()
    in_polished_writer.close()
    in_rejected_writer.close()
    out_polished_writer.close()

    def _f(bed_fn, what):
        n = len([r for r in BedReader(bed_fn)])
        log.info("There are {n} {what} structural variants in file {f} .".format(n=n, what=what, f=bed_fn))

    _f(in_passed_bed_fn, 'PASSED')
    _f(in_skipped_bed_fn, 'SKIPPED')
    _f(in_rejected_bed_fn, 'REJECTED')
    _f(in_polished_bed_fn, 'INPUT POLISHED')
    _f(out_polished_bed_fn, 'OUTPUT POLISHED')
    _f(collected_bed_fn, 'COLLECTED POLISHED')
