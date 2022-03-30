from .template import Processor
from .constant import TUMOR, NORMAL


class Vcf2Maf(Processor):

    SPECIES = 'homo_sapiens'
    NCBI_BUILD = 'GRCh38'

    annotated_vcf: str
    ref_fa: str

    maf: str

    def main(
            self,
            annotated_vcf: str,
            ref_fa: str) -> str:

        self.annotated_vcf = annotated_vcf
        self.ref_fa = ref_fa

        self.maf = f'{self.outdir}/variants.maf'
        cmd = self.CMD_LINEBREAK.join([
            'vcf2maf.pl',
            f'--input-vcf {self.annotated_vcf}',
            f'--ref-fasta {self.ref_fa}',
            f'--tmp-dir {self.workdir}',
            f'--output-maf {self.maf}',
            f'--tumor-id {TUMOR}',
            f'--normal-id {NORMAL}',
            f'--vcf-tumor-id {TUMOR}',
            f'--vcf-normal-id {NORMAL}',
            '--inhibit-vep',
            f'--species {self.SPECIES}',
            f'--ncbi-build {self.NCBI_BUILD}',
        ])
        self.call(cmd)

        return self.maf


