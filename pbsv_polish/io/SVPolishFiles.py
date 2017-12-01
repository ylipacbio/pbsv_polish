#!/usr/bin/env python
import os.path as op
from ..independent import Constants as C
from pbsv.independent.utils import execute
from ..utils import write_to_bash_file, qsub_to_sge_or_run_local
from ..independent.cmds import (sv_pbdagcon_cmd, blasr_cmd,
        sort_index_bam_inline_cmd, pbindex_cmd, variant_caller_cmd,
        trim_lq_cmd, pbsv_run_and_transform_cmds)


class SVPolishFiles(object):

    def __init__(self, root_dir, min_qv, ref_ext_len):
        """
        root_dir is a directory containing all files for polishing structural variants.
        min_qv: minimimum polish QV to accept a polished SV.
        ref_ext_len: extend reference by n base pairs.
        """
        def f(strs):
            return '.'.join(strs)

        def g(fn):
            return op.join(self.root_dir, fn)

        self.root_dir = root_dir
        self.min_qv = min_qv
        self.ref_ext_len = ref_ext_len

        self.subreads_consensus_sh = g('subreads_consensus.sh')
        self.readme = g('README')

        self._subreads_prefix = g('sr')
        self.subreads_bam = f([self._subreads_prefix, 'subreads.bam'])
        self.subreads_fa = f([self._subreads_prefix, 'subreads.fasta'])

        self._dagcon_prefix = g('sv_pbdagcon')
        self.dagcon_fa = f([self._dagcon_prefix, 'fasta'])

        self._sv_ref_prefix = g('sv_ref_w_ext')
        self.sv_ref_fa = f([self._sv_ref_prefix, 'fasta'])

        self.sr_dagcon_blasr_bam = g('sr.sv_pbdagcon.blasr.bam')

        self._polish_prefix = g('polish.qv%s' % self.min_qv)
        self.polish_fa = f([self._polish_prefix, 'fasta'])
        self.polish_fq = f([self._polish_prefix, 'fastq'])
        self.polish_blasr_bed = f([self._polish_prefix, 'blasr', 'bed'])
        self.polish_blasr_sh = f([self._polish_prefix, 'blasr', 'sh'])
        self.polish_ref_blasr_bam = f([self._polish_prefix, 'ref', 'blasr', 'bam'])

        self.polish_ngmlr_bed = f([self._polish_prefix, 'ngmlr', 'bed'])
        self.polish_ngmlr_sh = f([self._polish_prefix, 'ngmlr', 'sh'])
        self.polish_ref_ngmlr_bam = f([self._polish_prefix, 'ref', 'ngmlr', 'bam'])

        self._polish_hqlq_prefix = g('polish.hqlq')
        self.polish_hqlq_fa = f([self._polish_hqlq_prefix, 'fasta'])
        self.polish_hqlq_fq = f([self._polish_hqlq_prefix, 'fastq'])

        self.run_sh = g('run.sh')

    def make_readme(self):
        """Write all files to readme"""
        with open(self.readme, 'w') as writer:
            writer.write('\n'.join([self.polish_blasr_sh, self.polish_ngmlr_sh, self.subreads_consensus_sh]))

    @property
    def scripts(self):
        return [self.subreads_consensus_sh, self.polish_ngmlr_sh, self.polish_blasr_sh]

    def make_all_scripts(self):
        self.make_polish_ngmlr_script()
        self.make_polish_blasr_script()
        self.make_subreads_consensus_script()
        write_to_bash_file(cmds=self.scripts, bash_sh_fn=self.run_sh)
        for sh_fn in self.scripts + [self.run_sh]:
            execute('chmod +x %s' % sh_fn)

    def execute_all_scripts(self, use_sge):
        try:
            qsub_to_sge_or_run_local(script_fn=self.run_sh, use_sge=use_sge)
        except Exception:
            print 'FAILED TO POLISH SV %s' % self.root_dir

    def make_polish_ngmlr_script(self):
        write_to_bash_file(cmds=self.polish_ngmlr_cmds, bash_sh_fn=self.polish_ngmlr_sh)

    def make_polish_blasr_script(self):
        write_to_bash_file(cmds=self.polish_blasr_cmds, bash_sh_fn=self.polish_blasr_sh)

    def make_subreads_consensus_script(self):
        write_to_bash_file(cmds=self.subreads_consensus_cmds, bash_sh_fn=self.subreads_consensus_sh)

    @property
    def subreads_consensus_cmds(self):
        """Return a list of cmds used to make consensus sequence of subreads_bam.
        sv_prefix -- a prefix string from a BedRecord obj, e.g., chr1_0_100_Deletion_-100
        """
        align_bam = self.sr_dagcon_blasr_bam
        # c0 will generate sv_pbdagcon output consensus sequence of subreads with fai
        c0 = sv_pbdagcon_cmd(self.subreads_bam, self._dagcon_prefix, 'subreads_consensus', self.sv_ref_fa)
        c1 = blasr_cmd(query_fn=self.subreads_bam, target_fn=self.dagcon_fa, out_fn=align_bam, nproc=C.BLASR_NPROC)
        c2 = sort_index_bam_inline_cmd(align_bam)
        c3 = pbindex_cmd(align_bam)
        c4 = variant_caller_cmd(align_bam=align_bam, ref_fa=self.dagcon_fa,
                                out_fa=self.polish_hqlq_fa, out_fq=self.polish_hqlq_fq, nproc=C.VARIANT_CALLER_NPROC)
        c5 = trim_lq_cmd(in_fq=self.polish_hqlq_fq, min_qv=self.min_qv, out_fq=self.polish_fq, out_fa=self.polish_fa)
        return [c0, c1, c2, c3, c4, c5]

    @property
    def polish_ngmlr_cmds(self):
        """A list of shell commands to ngmlr align polished sequence to a substring of chromosome, call
        structural variants, and transform coordinate back to the original chromosome"""
        return pbsv_run_and_transform_cmds(reads_fn=self.polish_fa, ref_fa_fn=self.sv_ref_fa,
                                           cfg_fn=C.PBSV_POLISH_CFG, o_bam_fn=self.polish_ref_ngmlr_bam,
                                           o_bed_fn=self.polish_ngmlr_bed, algorithm='ngmlr')

    @property
    def polish_blasr_cmds(self):
        """A list of shell commands to blasr align polished sequence to a substring of chromosome, call
        structural variants, and transform coordinate back to the original chromosome"""
        return pbsv_run_and_transform_cmds(reads_fn=self.polish_fa, ref_fa_fn=self.sv_ref_fa,
                                           cfg_fn=C.PBSV_POLISH_CFG, o_bam_fn=self.polish_ref_blasr_bam,
                                           o_bed_fn=self.polish_blasr_bed, algorithm='blasr')
