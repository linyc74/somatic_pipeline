import os
from .template import Processor


class SnpEff(Processor):

    GENOME_ID = 'GRCh38.99'
    SUMMARY_PREFIX = 'snpeff-summary'

    vcf: str
    annotated_vcf: str

    def main(self, vcf: str) -> str:
        self.vcf = vcf
        self.execute()
        self.move_summary_files()
        return self.annotated_vcf

    def execute(self):
        self.annotated_vcf = f'{self.outdir}/annotated.vcf'
        html = f'{self.outdir}/{self.SUMMARY_PREFIX}.html'
        stderr = f'{self.outdir}/snpeff.log'

        cmd = self.CMD_LINEBREAK.join([
            'snpeff',
            '-verbose',
            f'-htmlStats {html}',
            self.GENOME_ID,
            self.vcf,
            f'> {self.annotated_vcf}',
            f'2> {stderr}',
        ])

        self.call(cmd)

    def move_summary_files(self):
        dstdir = f'{self.outdir}/snpeff'
        os.makedirs(dstdir, exist_ok=True)
        self.call(f'mv {self.outdir}/{self.SUMMARY_PREFIX}* {dstdir}/')
