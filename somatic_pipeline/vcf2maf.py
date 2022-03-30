from .template import Processor
from .constant import TUMOR, NORMAL


class Vcf2Maf(Processor):

    SPECIES = 'homo_sapiens'
    NCBI_BUILD = 'GRCh38'

    annotated_vcf: str
    ref_fa: str
    variant_caller: str

    vcf_tumor_column: str
    vcf_normal_column: str
    maf: str

    def main(
            self,
            annotated_vcf: str,
            ref_fa: str,
            variant_caller: str) -> str:

        self.annotated_vcf = annotated_vcf
        self.ref_fa = ref_fa
        self.variant_caller = variant_caller

        self.set_vcf_columns()
        self.execute()

        return self.maf

    def set_vcf_columns(self):
        if self.variant_caller == 'mutect2':
            # Mutect2 uses tagged read group names as vcf column names
            self.vcf_tumor_column = TUMOR
            self.vcf_normal_column = NORMAL
        else:
            # hard-coded by most variant callers, e.g. MuSE and VarScan
            self.vcf_tumor_column = 'TUMOR'
            self.vcf_normal_column = 'NORMAL'

    def execute(self):
        log = f'{self.outdir}/vcf2maf.log'
        self.maf = f'{self.outdir}/variants.maf'
        cmd = self.CMD_LINEBREAK.join([
            'vcf2maf.pl',
            f'--input-vcf {self.annotated_vcf}',
            f'--ref-fasta {self.ref_fa}',
            f'--tmp-dir {self.workdir}',
            f'--output-maf {self.maf}',
            f'--tumor-id {self.vcf_tumor_column}',
            f'--normal-id {self.vcf_normal_column}',
            f'--vcf-tumor-id {self.vcf_tumor_column}',
            f'--vcf-normal-id {self.vcf_normal_column}',
            '--inhibit-vep',
            f'--species {self.SPECIES}',
            f'--ncbi-build {self.NCBI_BUILD}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)
