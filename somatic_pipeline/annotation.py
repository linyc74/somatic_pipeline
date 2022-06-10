import os
from os.path import basename
from typing import Optional, List
from .tools import edit_fpath
from .template import Processor


class Annotation(Processor):

    SNPEFF = 'snpeff'
    VEP = 'vep'

    annotator: str
    vcf: str
    ref_fa: str
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]
    vep_db_tar_gz: Optional[str]
    vep_db_type: str
    cadd_resource: Optional[str]

    def main(
            self,
            annotator: str,
            vcf: str,
            ref_fa: str,
            clinvar_vcf_gz: Optional[str],
            dbsnp_vcf_gz: Optional[str],
            snpsift_dbnsfp_txt_gz: Optional[str],
            vep_db_tar_gz: Optional[str],
            vep_db_type: str,
            cadd_resource: Optional[str]) -> str:

        self.annotator = annotator
        self.vcf = vcf
        self.ref_fa = ref_fa
        self.clinvar_vcf_gz = clinvar_vcf_gz
        self.dbsnp_vcf_gz = dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = snpsift_dbnsfp_txt_gz
        self.vep_db_tar_gz = vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.cadd_resource = cadd_resource

        assert self.annotator in [self.SNPEFF, self.VEP]

        if self.annotator == self.SNPEFF:
            self.run_snpeff()
            self.run_snpsift_annotate()
            self.run_snpsift_dbnsfp()
        else:
            self.run_vep()

        self.move_vcf()

        return self.vcf

    def run_snpeff(self):
        self.vcf = SnpEff(self.settings).main(vcf=self.vcf)

    def run_snpsift_annotate(self):
        for resource in [self.clinvar_vcf_gz, self.dbsnp_vcf_gz]:
            if resource is not None:
                self.vcf = SnpSiftAnnotate(self.settings).main(vcf=self.vcf, db_vcf_gz=resource)

    def run_snpsift_dbnsfp(self):
        if self.snpsift_dbnsfp_txt_gz is not None:
            self.vcf = SnpSiftDbNSFP(self.settings).main(
                vcf=self.vcf,
                dbnsfp_txt_gz=self.snpsift_dbnsfp_txt_gz)

    def run_vep(self):
        assert self.vep_db_tar_gz is not None
        self.vcf = VEP(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.ref_fa,
            vep_db_tar_gz=self.vep_db_tar_gz,
            vep_db_type=self.vep_db_type,
            cadd_resource=self.cadd_resource)

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
        stderr = f'{self.outdir}/snpsift-annotate.log'
        cmd = self.CMD_LINEBREAK.join([
            'snpsift annotate',
            self.db_vcf_gz,
            self.vcf,
            f'> {self.output_vcf}',
            f'2>> {stderr}',
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
        stderr = f'{self.outdir}/snpsift-dbnsfp.log'
        cmd = self.CMD_LINEBREAK.join([
            'snpsift dbnsfp',
            f'-db {self.dbnsfp_txt_gz}',
            self.vcf,
            f'> {self.output_vcf}',
            f'2> {stderr}',
        ])
        self.call(cmd)


class VEP(Processor):

    VEP_DB_TYPE_TO_FLAG = {
        'merged': '--merged',
        'refseq': '--refseq',
        'vep': '',
    }

    vcf: str
    ref_fa: str
    vep_db_tar_gz: str
    vep_db_type: str
    cadd_resource: Optional[str]

    cache_dir: str
    plugin_args: List[str]
    output_vcf: str

    def main(
            self,
            vcf: str,
            ref_fa: str,
            vep_db_tar_gz: str,
            vep_db_type: str,
            cadd_resource: Optional[str]) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa
        self.vep_db_tar_gz = vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.cadd_resource = cadd_resource

        self.extract_db_tar_gz()
        self.set_plugin_args()
        self.set_output_vcf()
        self.execute()
        self.move_summary_html()

        return self.output_vcf

    def extract_db_tar_gz(self):
        self.cache_dir = f'{self.workdir}/.vep'
        os.makedirs(self.cache_dir, exist_ok=True)
        self.call(f'tar -xzf {self.vep_db_tar_gz} -C "{self.cache_dir}"')

    def set_plugin_args(self):
        self.plugin_args = []
        if self.cadd_resource is not None:
            self.plugin_args.append(f'--plugin CADD,{self.cadd_resource}')

    def set_output_vcf(self):
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-vep.vcf',
            dstdir=self.workdir)

    def execute(self):
        log = f'{self.outdir}/vep.log'
        args = [
            'vep',
            '--offline',
            f'--input_file {self.vcf}',
            f'--fasta {self.ref_fa}',
            f'--dir_cache {self.cache_dir}',
            self.VEP_DB_TYPE_TO_FLAG[self.vep_db_type],
        ] + self.plugin_args + [
            '--vcf',  # output in vcf format
            f'--output_file {self.output_vcf}',
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))

    def move_summary_html(self):
        src = f'{self.output_vcf}_summary.html'
        dst = f'{self.outdir}/vep-summary.html'
        self.call(f'mv {src} {dst}')
