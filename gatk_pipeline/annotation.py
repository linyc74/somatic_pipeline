from .constant import CMD_LINEBREAK
from .template import Processor, Settings


class SnpEff(Processor):

    GENOME_ID = 'GRCh38.99'

    vcf: str
    annotated_vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, vcf: str) -> str:
        self.vcf = vcf

        self.execute()

        return self.annotated_vcf

    def execute(self):
        self.annotated_vcf = f'{self.outdir}/annotated.vcf'
        html = f'{self.outdir}/snpEff_summary.html'
        stderr = f'{self.workdir}/snpEff_stderr.log'

        cmd = CMD_LINEBREAK.join([
            'snpeff',
            '-verbose',
            f'-htmlStats {html}',
            self.GENOME_ID,
            self.vcf,
            f'1> {self.annotated_vcf}',
            f'2> {stderr}',
        ])

        self.call(cmd)
