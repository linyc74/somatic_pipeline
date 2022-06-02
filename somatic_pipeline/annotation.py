import os
from os.path import basename
from typing import Optional
from .tools import edit_fpath
from .template import Processor


class Annotation(Processor):

    vcf: str
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]

    def main(
            self,
            vcf: str,
            clinvar_vcf_gz: Optional[str],
            dbsnp_vcf_gz: Optional[str],
            snpsift_dbnsfp_txt_gz: Optional[str]) -> str:

        self.vcf = vcf
        self.clinvar_vcf_gz = clinvar_vcf_gz
        self.dbsnp_vcf_gz = dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = snpsift_dbnsfp_txt_gz

        self.run_snpeff()
        self.run_snpsift_annotate()
        self.run_snpsift_dbnsfp()
        self.move_vcf()

        return self.vcf

    def run_snpeff(self):
        self.vcf = SnpEff(self.settings).main(vcf=self.vcf)

    def run_snpsift_annotate(self):
        for resource in [self.clinvar_vcf_gz, self.dbsnp_vcf_gz]:
            if resource is not None:
                SnpSiftAnnotate(self.settings).main(vcf=self.vcf, db_vcf_gz=resource)

    def run_snpsift_dbnsfp(self):
        if self.snpsift_dbnsfp_txt_gz is not None:
            self.vcf = SnpSiftDbNSFP(self.settings).main(
                vcf=self.vcf,
                dbnsfp_txt_gz=self.snpsift_dbnsfp_txt_gz)

    def move_vcf(self):
        dst = f'{self.outdir}/annotated.vcf'
        self.call(f'mv {self.vcf} {dst}')
        self.vcf = dst


class SnpEff(Processor):

    GENOME_ID = 'GRCh38.99'
    SUMMARY_PREFIX = 'snpeff-summary'

    vcf: str
    output_vcf: str

    def main(self, vcf: str) -> str:
        self.vcf = vcf
        self.set_output_vcf()
        self.execute()
        self.move_summary_files()
        return self.output_vcf

    def set_output_vcf(self):
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-snpeff.vcf',
            dstdir=self.workdir)

    def execute(self):
        html = f'{self.outdir}/{self.SUMMARY_PREFIX}.html'
        stderr = f'{self.outdir}/snpeff.log'
        cmd = self.CMD_LINEBREAK.join([
            'snpeff',
            '-verbose',
            f'-htmlStats {html}',
            self.GENOME_ID,
            self.vcf,
            f'> {self.output_vcf}',
            f'2> {stderr}',
        ])
        self.call(cmd)

    def move_summary_files(self):
        dstdir = f'{self.outdir}/snpeff'
        os.makedirs(dstdir, exist_ok=True)
        self.call(f'mv {self.outdir}/{self.SUMMARY_PREFIX}* {dstdir}/')


class SnpSiftAnnotate(Processor):

    vcf: str
    db_vcf_gz: str

    output_vcf: str

    def main(self, vcf: str, db_vcf_gz: str) -> str:

        self.vcf = vcf
        self.db_vcf_gz = db_vcf_gz

        self.prepare_db_vcf_gz()
        self.set_output_vcf()
        self.execute()

        return self.output_vcf

    def prepare_db_vcf_gz(self):
        self.__copy()
        self.__index()

    def __copy(self):
        src = self.db_vcf_gz
        dst = edit_fpath(
            fpath=self.db_vcf_gz,
            dstdir=self.workdir)
        self.call(f'cp {src} {dst}')
        self.db_vcf_gz = dst

    def __index(self):
        self.call(f'tabix -p vcf {self.db_vcf_gz}')

    def set_output_vcf(self):
        suffix = basename(self.db_vcf_gz).rstrip('.vcf.gz')
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix=f'-{suffix}.vcf',
            dstdir=self.workdir)

    def execute(self):
        stderr = f'{self.outdir}/snpsift.log'
        cmd = self.CMD_LINEBREAK.join([
            'snpsift annotate',
            self.db_vcf_gz,
            self.vcf,
            f'> {self.output_vcf}',
            f'2> {stderr}',
        ])
        self.call(cmd)


class SnpSiftDbNSFP(Processor):

    vcf: str
    dbnsfp_txt_gz: str

    output_vcf: str

    def main(
            self,
            vcf: str,
            dbnsfp_txt_gz: str) -> str:
        self.vcf = vcf
        self.dbnsfp_txt_gz = dbnsfp_txt_gz

        self.prepare_dbnsfp_txt_gz()
        self.set_output_vcf()
        self.execute()

        return self.output_vcf

    def prepare_dbnsfp_txt_gz(self):
        self.__copy()
        self.__index()

    def __copy(self):
        src = self.dbnsfp_txt_gz
        dst = edit_fpath(
            fpath=self.dbnsfp_txt_gz,
            dstdir=self.workdir)
        self.call(f'cp {src} {dst}')
        self.dbnsfp_txt_gz = dst

    def __index(self):
        self.call(f'tabix -p vcf {self.dbnsfp_txt_gz}')

    def set_output_vcf(self):
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-dbnsfp.vcf',
            dstdir=self.workdir)

    def execute(self):
        stderr = f'{self.outdir}/snpsift.log'
        cmd = self.CMD_LINEBREAK.join([
            'snpsift dbnsfp',
            f'-db {self.dbnsfp_txt_gz}',
            self.vcf,
            f'> {self.output_vcf}',
            f'2> {stderr}',
        ])
        self.call(cmd)
