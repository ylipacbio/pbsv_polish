#!/usr/bin/env python
import os.path as op
import sys
import os
from collections import defaultdict
from pbcore.io import DataSet, FastaWriter

from pbsv.independent.utils import execute, realpath, write_to_bash_file, execute_as_bash

from .utils import *


def make_consensus_script_of_subreads_bam(subreads_bam, o_script_fn, o_consensus_id, o_consensus_fn):
    """Make a script file which when excuted, make a consensus sequence of subreads_bam"""
    cmds = []
    output_prefix = op.join(op.dirname(o_consensus_fn), 'sv_pbdagcon')
    c0 = 'sv_pbdagcon %s %s %s' % (subreads_bam, output_prefix, o_consensus_id)
    output_dagcon_fasta = output_prefix + '_ref.fasta'
    align_bam = op.join(op.dirname(o_consensus_fn), 'sr2_sv_pbdagcon.bam')
    nproc = 16
    c1 = 'blasr %s %s --nproc 16 --hitPolicy randomBest --minMatch 8 --maxMatch 15 --bam --out %s' % (subreads_bam,output_dagcon_fasta, align_bam)
    c2 = 'samtools faidx %s' % (output_dagcon_fasta)
    c3 = 'samtools sort {f1} -o {f2} && mv {f2} {f1} && samtools index {f1} && pbindex {f1}'.format(f1=align_bam, f2=align_bam+'.sorted.bam')
    output_prefix = o_consensus_fn[0:o_consensus_fn.rfind('.')]
    hqlq_out_fa, hqlq_out_fq = output_prefix + '.hqlq.fasta', output_prefix + '.hqlq.fastq'
    out_fa, out_fq = output_prefix + '.fasta', output_prefix + '.fastq'
    c4 = "variantCaller --algorithm=best {aln_bam} --verbose -j{nproc} --minMapQV {minqv} --referenceFilename={ref_fa} -o {out_fa} -o {out_fq}".format(aln_bam=align_bam, ref_fa=output_dagcon_fasta, out_fa=hqlq_out_fa, out_fq=hqlq_out_fq, nproc=nproc, minqv=10)
    c5 = 'trim_lq %s %s --min_qv 20' % (hqlq_out_fq, out_fq) # simply remove lower case sequences on both ends
    c6 = 'fq2fa %s %s' % (out_fq, out_fa)
    cmds = [c0, c1, c2, c3, c4, c5, c6]

    print 'Running %s' % o_script_fn
    write_to_bash_file(cmds, o_script_fn)
    return o_script_fn

def qsub_to_sge_or_run_local(script_fn, use_sge=False):
    if use_sge:
        cmd = 'qsub -q def66 -pe smp 1 -v -cwd -b y {f} 1>{f}.stdout.log 2>{f}.stderr.log'.format(f=script_fn)
    else:
        cmd = 'chmod +x %s && bash %s' % (script_fn, script_fn)
    execute(cmd)

def make_script_of_pbsv_ngmlr_call(reads_fn, ref_fasta_fn, cfg_fn, o_prefix):
    cmds = []
    o_sv_fn = o_prefix + '.bed'
    o_script_fn = o_prefix + '.sh'
    tmp_sv_fn  = o_script_fn + '.use_substr_as_chrom.bed'
    chained_bam_fn = o_prefix + '.chained.bam'
    c0 = 'pbsv align %s %s %s --cfg_fn %s' % (ref_fasta_fn, reads_fn, chained_bam_fn, cfg_fn)
    c1 = 'pbsv call %s %s %s --cfg_fn %s' % (ref_fasta_fn, chained_bam_fn, tmp_sv_fn, cfg_fn)
    c2 = 'sv_transform_coordinate %s %s' % (tmp_sv_fn, o_sv_fn)
    cmds = [c0, c1, c2]
    print 'Running %s' % o_script_fn
    write_to_bash_file(cmds, o_script_fn)
    return o_script_fn

def make_script_of_pbsv_blasr_call(reads_fn, ref_fasta_fn, cfg_fn, o_prefix):
    """Using blasr to make alignments."""
    cmds = []
    o_sv_fn = o_prefix + '.bed'
    o_script_fn = o_prefix + '.sh'
    tmp_sv_fn  = o_script_fn + '.use_substr_as_chrom.bed'
    chained_bam_fn = o_prefix + '.chained.bam'
    c0 = 'blasr %s %s --bam --out %s --bestn 1 --maxMatch 15' % (reads_fn, ref_fasta_fn, chained_bam_fn)
    c1 = 'pbsv call %s %s %s --cfg_fn %s' % (ref_fasta_fn, chained_bam_fn, tmp_sv_fn, cfg_fn)
    c2 = 'sv_transform_coordinate %s %s' % (tmp_sv_fn, o_sv_fn)
    c3 = 'samtools sort {f} -o {f}.tmp && mv {f}.tmp {f} && samtools index {f}'.format(f=chained_bam_fn)
    cmds = [c0, c1, c2, c3]
    print 'Running %s' % o_script_fn
    write_to_bash_file(cmds, o_script_fn)
    return o_script_fn

def make_diagnose_script_for_pbsv_run(o_dir):
    diagnose_fn = op.join(o_dir, 'diagnose.sh')
    cmd = """
fastalen polished.fasta
fastalen polished.hqlq.fasta
blasr polished.fasta sv_reference_w_extension.fasta --header --maxMatch 15 -m 4
blasr polished.hqlq.fasta sv_reference_w_extension.fasta --header --maxMatch 15 -m 4
"""
    with open(diagnose_fn, 'w') as w:
        w.write(cmd)

def write_fasta(o_fasta_fn, records):
    """Write a list of fasta records [(name, seq), ...,  (name, seq)] to o_fasta_fn"""
    with FastaWriter(o_fasta_fn) as w:
        for r in records:
            w.writeRecord(r[0], r[1])

def substr_fasta(fileobj, chrom, start, end, o_fasta_fn):
    """fetch a substring of reference fasta sequence and save to output fasta file o_fasta_fn"""
    try:
        seq = fileobj.fetch(chrom, start, end)
    except Exception as e:
        raise ValueError("Could not get substring (%s, %s, %s) from %s" % (chrom, start, end, fileobj.filename))
    name = '%s__substr__%s_%s' % (chrom, start, end)
    write_fasta(o_fasta_fn, [(name, seq)])


def main():
    #in_dir = "in_sl_228_hg00733_10fold"
    in_dir = 'in_sl_1813_yeast_10fold'
    aln_fn = op.join(in_dir, "alignments.bam")
    subreads_xml_fn = op.join(in_dir, "subreads.xml")
    genome_fa = op.join(in_dir, "genome.fa")
    bed_fn = op.join(in_dir, 'structural_variants.bed')
    #bed_fn = op.join(in_dir, 'chrV_116286_116286_Insertion_5899.bed')
    out_dir = 'out'
    reference_fasta_obj = Fastafile(genome_fa)
    REFERENCE_EXTENSION = 200000

    ofile_obj = open('coverage.txt', 'w')
    alnfile_obj = X2PysamReader(aln_fn)._alignment_file
    bedreader_obj = BedReader(bed_fn)
    movie2bams = get_movies2bams_from_subreads_xml(subreads_xml_fn)
    print movie2bams

    POLISH_CFG_FN = '/pbi/dept/secondary/siv/yli/sv/consensus/pbsv.polish.cfg'
    make_bam = False
    i = 0
    for bed_record, alns in yield_alns_from_bed_file(alnfile_obj, bedreader_obj=bedreader_obj):
        i += 1
        if i % 1000 == 0:
            print "i=%s, got %s covering alignemnts for sv %s" % (i, len(alns), bed_record)
        srs = get_query_subreads_from_alns(alns)
        zmws = get_query_zmws_from_alns(alns)
        ofile_obj.write('%s\t%s\t%s\n' % (len(srs), len(zmws), bed_record))

        if True:
            # make a subdirectory
            sv_str = '_'.join([str(x) for x in [bed_record.chrom, bed_record.start, bed_record.end, bed_record.sv_type, bed_record.sv_len]])
            data_dir = realpath(op.join(out_dir, sv_str))
            _mkdir(data_dir)
            out_prefix = op.join(data_dir, 'in')
            subreads_bam = make_subreads_bam_of_zmws(movie2bams=movie2bams, zmws=zmws, out_prefix=out_prefix, dry_run=(not make_bam))
            script_fn = op.join(data_dir, 'make_consensus.sh')
            consensus_fn = op.join(data_dir, 'polished.fasta')
            f0 = make_consensus_script_of_subreads_bam(subreads_bam=subreads_bam, o_script_fn=script_fn, o_consensus_id=sv_str, o_consensus_fn=consensus_fn)

            ref_start, ref_end = max(0, bed_record.start - REFERENCE_EXTENSION), bed_record.end + REFERENCE_EXTENSION
            ref_fasta_fn = op.join(data_dir, 'sv_reference_w_extension.fasta')
            substr_fasta(fileobj=reference_fasta_obj, chrom=bed_record.chrom, start=ref_start, end=ref_end, o_fasta_fn=ref_fasta_fn)

            o_prefix = op.join(data_dir, 'polish.qv20.ngmlr')
            f1 = make_script_of_pbsv_ngmlr_call(reads_fn=consensus_fn, ref_fasta_fn=ref_fasta_fn, cfg_fn=POLISH_CFG_FN, o_prefix=o_prefix)

            o_prefix = op.join(data_dir, 'polish.hqlq.ngmlr')
            consensus_hqlq_fa = op.join(data_dir, 'polished.hqlq.fasta')
            f2 = make_script_of_pbsv_ngmlr_call(reads_fn=consensus_hqlq_fa, ref_fasta_fn=ref_fasta_fn, cfg_fn=POLISH_CFG_FN, o_prefix=o_prefix)

            o_prefix = op.join(data_dir, 'polish.qv20.blasr')
            f3 = make_script_of_pbsv_blasr_call(reads_fn=consensus_fn, ref_fasta_fn=ref_fasta_fn, cfg_fn=POLISH_CFG_FN, o_prefix=o_prefix)

            #o_prefix = op.join(data_dir, 'polish.hqlq.blasr')
            #make_script_of_pbsv_blasr_call(reads_fn=consensus_hqlq_fa, ref_fasta_fn=ref_fasta_fn, cfg_fn=POLISH_CFG_FN, o_prefix=o_prefix)

        #if True:
        #    sv_str = '_'.join([str(x) for x in [bed_record.chrom, bed_record.start, bed_record.end, bed_record.sv_type, bed_record.sv_len]])
        #    data_dir = realpath(op.join(out_dir, sv_str))
        #    f1 = op.join(data_dir, 'polish.qv20.ngmlr.sh')
        #    f2 = op.join(data_dir, 'polish.hqlq.ngmlr.sh')
        #    f3 = op.join(data_dir, 'polish.qv20.blasr.sh')

            sh_fns = [f1, f2, f3] # don't run f0 consensus to save time debugging
            f_all = op.join(data_dir, 'all_script.sh')
            cmds = ['bash %s' % f for f in sh_fns]
            print 'qsub %s' % f_all
            write_to_bash_file(cmds, f_all)
            qsub_to_sge_or_run_local(f_all, use_sge=True)

            make_diagnose_script_for_pbsv_run(o_dir=data_dir)
            if i == 1000:
                break

    ofile_obj.close()

if __name__ == "__main__":
    main()
