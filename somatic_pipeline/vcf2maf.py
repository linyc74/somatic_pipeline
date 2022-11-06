from os.path import dirname
from typing import Optional
from .tools import edit_fpath
from .template import Processor


class Vcf2Maf(Processor):

    SPECIES = 'homo_sapiens'
    NCBI_BUILD = 'GRCh38'

    vcf: str
    ref_fa: str
    dstdir: Optional[str]

    vcf_tumor_column: str
    vcf_normal_column: str
    maf: str

    def main(
            self,
            vcf: str,
            ref_fa: str,
            dstdir: Optional[str]) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa
        self.dstdir = dstdir

        self.set_tumor_normal_columns()
        self.set_maf()
        self.execute()

        return self.maf

    def set_tumor_normal_columns(self):
        # Default values
        self.vcf_tumor_column = 'T'
        self.vcf_normal_column = 'N'

        with open(self.vcf) as fh:
            for line in fh:
                if line.startswith('#CHROM'):
                    last_two_columns = line.strip().split('\t')[-2:]
                    break

        for c in last_two_columns:
            if c.lower() == 'tumor':
                self.vcf_tumor_column = c
            if c.lower() == 'normal':
                self.vcf_normal_column = c

    def set_maf(self):
        self.maf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='.maf',
            dstdir=dirname(self.vcf) if self.dstdir is None else self.dstdir
        )

    def execute(self):
        log = f'{self.outdir}/vcf2maf.log'
        cmd = self.CMD_LINEBREAK.join([
            'vcf2maf.pl',
            f'--input-vcf {self.vcf}',
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
