import pytest

import os.path as op
from pbsv_polish.cli import  pbsvp_main, pbsvputil_main
from setpath import SIV_DATA_DIR, OUT_DIR
from pbsv.independent.utils import execute, rmpath


def test_pbsvp(capsys):
    # Test pbsvp polish
    in_data_dir = op.join(SIV_DATA_DIR, 'FY1679-12fold')
    genome_fa, subreads_xml, alignments_bam, rich_bed = [op.join(in_data_dir, f) for f in
            ['genome.fa', 'subreads.xml', 'alignments.bam', 'sv.rich.bed']]
    work_dir = op.join(OUT_DIR, 'test-pbsvp')

    rmpath(work_dir)
    pbsvp_main( ['main', 'polish', genome_fa, subreads_xml, alignments_bam, rich_bed, work_dir])

    for f in ['in.passed.bed', 'in.skipped.bed', 'chr13_202206_202206_Insertion_152']:
        assert op.exists(op.join(work_dir, f))

    # Execute
    execute('bash {}'.format(op.join(work_dir, 'chr13_202206_202206_Insertion_152', 'run.sh')))

    # Test pbsvp collect
    pbsvp_main(['main', 'collect', work_dir, op.join(work_dir, 'collected.bed')])

    for f in ['in.rejected.bed', 'in.polished.bed', 'out.polished.bed']:
        assert op.exists(op.join(work_dir, f))
