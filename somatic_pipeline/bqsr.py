from .tools import edit_fpath
from .template import Processor
from .index_files import SamtoolsIndexFa, GATKCreateSequenceDictionary, GATKIndexVcf


class BQSR(Processor):

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

        self.prepare_ref_fa()
        self.prepare_known_variant_vcf()
        self.base_recalibrator()
        self.apply_bqsr()

        return self.out_bam

    def prepare_ref_fa(self):
        SamtoolsIndexFa(self.settings).main(fa=self.ref_fa)
        GATKCreateSequenceDictionary(self.settings).main(ref_fa=self.ref_fa)

    def prepare_known_variant_vcf(self):
        src = self.known_variant_vcf
        dst = edit_fpath(
            fpath=self.known_variant_vcf,
            dstdir=self.workdir)
        self.call(f'cp {src} {dst}')
        self.known_variant_vcf = dst

        GATKIndexVcf(self.settings).main(vcf=self.known_variant_vcf)

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
