from .template import Processor


class SnpEff(Processor):

    GENOME_ID = 'GRCh38.99'

    vcf: str
    annotated_vcf: str

    def main(self, vcf: str) -> str:
        self.vcf = vcf

        self.execute()

        return self.annotated_vcf

    def execute(self):
        self.annotated_vcf = f'{self.outdir}/annotated.vcf'
        html = f'{self.outdir}/snpEff_summary.html'
        stderr = f'{self.outdir}/snpEff.log'

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
