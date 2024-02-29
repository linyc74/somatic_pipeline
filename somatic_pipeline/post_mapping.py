import os
from os.path import samefile, dirname
from typing import Optional, Tuple
from .tools import edit_fpath
from .template import Processor
from .index_files import SamtoolsIndexFa, GATKCreateSequenceDictionary, GATKIndexVcf


class MarkDuplicates(Processor):

    tumor_bam: str
    normal_bam: Optional[str]

    out_tumor_bam: str
    out_normal_bam: Optional[str]

    def main(
            self,
            tumor_bam: str,
            normal_bam: Optional[str]) -> Tuple[str, Optional[str]]:

        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        self.out_tumor_bam = GATKMarkDuplicates(self.settings).main(
            bam=self.tumor_bam)

        self.out_normal_bam = None if self.normal_bam is None else \
            GATKMarkDuplicates(self.settings).main(bam=self.normal_bam)

        return self.out_tumor_bam, self.out_normal_bam


class GATKMarkDuplicates(Processor):

    REMOVE_DUPLICATES = 'false'
    METRICS_DIRNAME = 'duplicate-metrics'

    bam: str

    metrics_txt: str
    out_bam: str

    def main(self, bam: str) -> str:
        self.bam = bam
        self.set_out_bam()
        self.set_metrics_txt()
        self.execute()
        return self.out_bam

    def set_out_bam(self):
        self.out_bam = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='-mark-duplicates.bam',
            dstdir=self.workdir)

    def set_metrics_txt(self):
        dstdir = f'{self.outdir}/{self.METRICS_DIRNAME}'
        os.makedirs(dstdir, exist_ok=True)
        self.metrics_txt = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='-duplicate-metrics.txt',
            dstdir=dstdir)

    def execute(self):
        log = f'{self.outdir}/gatk-MarkDuplicates.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk MarkDuplicates',
            f'--INPUT {self.bam}',
            f'--METRICS_FILE {self.metrics_txt}',
            f'--OUTPUT {self.out_bam}',
            f'--REMOVE_DUPLICATES {self.REMOVE_DUPLICATES}',
            f'1> {log} 2> {log}',
        ])
        self.call(cmd)


#


class BQSR(Processor):

    tumor_bam: str
    normal_bam: Optional[str]
    ref_fa: str
    known_variant_vcf: str

    temp_known_variant_vcf: str

    def main(
            self,
            tumor_bam: str,
            normal_bam: Optional[str],
            ref_fa: str,
            known_variant_vcf: str) -> Tuple[str, str]:

        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.ref_fa = ref_fa
        self.known_variant_vcf = known_variant_vcf

        self.prepare_ref_fa()
        self.prepare_known_variant_vcf()

        self.tumor_bam = self.run_bqsr(self.tumor_bam)
        if self.normal_bam is not None:
            self.normal_bam = self.run_bqsr(self.normal_bam)

        self.call(f'rm {self.temp_known_variant_vcf}')  # remove temp file to save disk space

        return self.tumor_bam, self.normal_bam

    def prepare_ref_fa(self):
        SamtoolsIndexFa(self.settings).main(fa=self.ref_fa)
        GATKCreateSequenceDictionary(self.settings).main(ref_fa=self.ref_fa)

    def prepare_known_variant_vcf(self):
        self.temp_known_variant_vcf = edit_fpath(
            fpath=self.known_variant_vcf,
            dstdir=self.workdir)

        self.call(f'cp {self.known_variant_vcf} {self.temp_known_variant_vcf}')

        GATKIndexVcf(self.settings).main(vcf=self.temp_known_variant_vcf)

    def run_bqsr(self, bam: str) -> str:
        return RunBQSR(self.settings).main(
            bam=bam,
            ref_fa=self.ref_fa,
            known_variant_vcf=self.temp_known_variant_vcf)


class RunBQSR(Processor):

    bam: str
    ref_fa: str
    known_variant_vcf: str

    recalibration_table: str
    out_bam: str

    def main(
            self,
            bam: str,
            ref_fa: str,
            known_variant_vcf: str) -> str:

        self.bam = bam
        self.ref_fa = ref_fa
        self.known_variant_vcf = known_variant_vcf

        self.base_recalibrator()
        self.apply_bqsr()
        self.remove_input_bam()  # the bam file before BQSR is no longer needed

        return self.out_bam

    def base_recalibrator(self):
        self.recalibration_table = f'{self.workdir}/BQSR-recalibration.table'
        log = f'{self.outdir}/gatk-BaseRecalibrator.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk BaseRecalibrator',
            f'--input {self.bam}',
            f'--reference {self.ref_fa}',
            f'--known-sites {self.known_variant_vcf}',
            f'--output {self.recalibration_table}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)

    def apply_bqsr(self):
        self.out_bam = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='-bqsr.bam',
            dstdir=self.workdir)
        log = f'{self.outdir}/gatk-ApplyBQSR.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk ApplyBQSR',
            f'--input {self.bam}',
            f'--reference {self.ref_fa}',
            f'--bqsr-recal-file {self.recalibration_table}',
            f'--output {self.out_bam}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)

    def remove_input_bam(self):
        if samefile(dirname(self.bam), self.workdir):  # if the input bam is in the workdir
            self.call(f'rm {self.bam}')


#


class MappingStats(Processor):

    tumor_bam: str
    normal_bam: Optional[str]

    def main(
            self,
            tumor_bam: str,
            normal_bam: Optional[str]):

        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        SamtoolsStats(self.settings).main(self.tumor_bam)
        if self.normal_bam is not None:
            SamtoolsStats(self.settings).main(self.normal_bam)


class SamtoolsStats(Processor):

    DSTDIR_NAME = 'mapping-stats'

    bam: str

    dstdir: str

    def main(self, bam: str):
        self.bam = bam
        self.make_dstdir()
        self.execute()

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def execute(self):
        txt = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='.txt',
            dstdir=self.dstdir)
        self.call(f'samtools stats {self.bam} > {txt}')
