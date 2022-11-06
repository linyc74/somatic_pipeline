import os
from os.path import basename
from typing import Optional, List
from .tools import edit_fpath
from .template import Processor


class Annotation(Processor):

    VEP = 'vep'
    SNPEFF = 'snpeff'

    annotator: str
    vcf: str
    ref_fa: str
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]
    vep_db_tar_gz: Optional[str]
    vep_db_type: str
    vep_buffer_size: int
    cadd_resource: Optional[str]
    dbnsfp_resource: Optional[str]

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
            vep_buffer_size: int,
            cadd_resource: Optional[str],
            dbnsfp_resource: Optional[str]) -> str:

        self.annotator = annotator
        self.vcf = vcf
        self.ref_fa = ref_fa
        self.clinvar_vcf_gz = clinvar_vcf_gz
        self.dbsnp_vcf_gz = dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = snpsift_dbnsfp_txt_gz
        self.vep_db_tar_gz = vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.vep_buffer_size = vep_buffer_size
        self.cadd_resource = cadd_resource
        self.dbnsfp_resource = dbnsfp_resource

        assert self.annotator in [self.SNPEFF, self.VEP]

        if self.annotator == self.SNPEFF:
            self.run_snpeff()
            self.run_snpsift_dbnsfp()
        else:
            self.run_vep()

        self.annotate_by_vcf_gz()  # independent of the choice of annotator
        self.move_vcf_to_outdir()

        return self.vcf

    def run_snpeff(self):
        self.vcf = SnpEff(self.settings).main(vcf=self.vcf)

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
            vep_buffer_size=self.vep_buffer_size,
            cadd_resource=self.cadd_resource,
            dbnsfp_resource=self.dbnsfp_resource)

    def annotate_by_vcf_gz(self):
        for vcf_gz in [self.clinvar_vcf_gz, self.dbsnp_vcf_gz]:
            if vcf_gz is not None:
                self.vcf = SnpSiftAnnotate(self.settings).main(vcf=self.vcf, resource_vcf_gz=vcf_gz)

    def move_vcf_to_outdir(self):
        dst = f'{self.outdir}/variants.vcf'
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
    resource_vcf_gz: str

    output_vcf: str

    def main(self, vcf: str, resource_vcf_gz: str) -> str:

        self.vcf = vcf
        self.resource_vcf_gz = resource_vcf_gz

        self.prepare_resource()
        self.set_output_vcf()
        self.execute()

        return self.output_vcf

    def prepare_resource(self):
        self.resource_vcf_gz = CopyAndTabixVcfGz(self.settings).main(self.resource_vcf_gz)

    def set_output_vcf(self):
        suffix = basename(self.resource_vcf_gz).rstrip('.vcf.gz')
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix=f'-{suffix}.vcf',
            dstdir=self.workdir)

    def execute(self):
        stderr = f'{self.outdir}/snpsift-annotate.log'
        cmd = self.CMD_LINEBREAK.join([
            'snpsift annotate',
            self.resource_vcf_gz,
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
        self.dbnsfp_txt_gz = CopyAndTabixTsvGz(self.settings).main(
            file=self.dbnsfp_txt_gz,
            seqname_column=1,
            start_column=2,
            end_column=2)

    def set_output_vcf(self):
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-snpsift-dbnsfp.vcf',
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
    NUM_THREADS = 4  # recommended on the VEP website

    vcf: str
    ref_fa: str
    vep_db_tar_gz: str
    vep_db_type: str
    vep_buffer_size: int
    cadd_resource: Optional[str]
    dbnsfp_resource: Optional[str]

    cache_dir: str
    plugin_args: List[str]
    output_vcf: str

    def main(
            self,
            vcf: str,
            ref_fa: str,
            vep_db_tar_gz: str,
            vep_db_type: str,
            vep_buffer_size: int,
            cadd_resource: Optional[str],
            dbnsfp_resource: Optional[str]) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa
        self.vep_db_tar_gz = vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.vep_buffer_size = vep_buffer_size
        self.cadd_resource = cadd_resource
        self.dbnsfp_resource = dbnsfp_resource

        self.extract_db_tar_gz()
        self.prepare_cadd_resource()
        self.prepare_dbnsfp_resource()
        self.set_plugin_args()
        self.set_output_vcf()
        self.execute()
        self.move_summary_html()

        return self.output_vcf

    def extract_db_tar_gz(self):
        self.cache_dir = f'{self.workdir}/.vep'
        os.makedirs(self.cache_dir, exist_ok=True)
        self.call(f'tar -xzf {self.vep_db_tar_gz} -C "{self.cache_dir}"')

    def prepare_cadd_resource(self):
        if self.cadd_resource is not None:
            self.cadd_resource = CopyAndTabixTsvGz(self.settings).main(
                file=self.cadd_resource,
                seqname_column=1,
                start_column=2,
                end_column=2)

    def prepare_dbnsfp_resource(self):
        if self.dbnsfp_resource is not None:
            AssertDbnsfpResourceFilenameForVEP(self.settings).main(self.dbnsfp_resource)
            self.dbnsfp_resource = CopyAndTabixTsvGz(self.settings).main(
                file=self.dbnsfp_resource,
                seqname_column=1,
                start_column=2,
                end_column=2)

    def set_plugin_args(self):
        self.plugin_args = []
        if self.cadd_resource is not None:
            self.plugin_args.append(f'--plugin CADD,{self.cadd_resource}')
        if self.dbnsfp_resource is not None:
            self.plugin_args.append(f'--plugin dbNSFP,{self.dbnsfp_resource},ALL')

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
            '--everything',
            f'--fork {self.threads}',
            f'--buffer_size {self.vep_buffer_size}',
        ]
        args += self.plugin_args
        args += [
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


class CopyAndTabix(Processor):

    file: str
    output: str

    def copy(self):
        self.output = edit_fpath(
            fpath=self.file,
            dstdir=self.workdir)
        self.call(f'cp {self.file} {self.output}')


class CopyAndTabixVcfGz(CopyAndTabix):

    def main(self, file: str) -> str:
        self.file = file
        self.copy()
        self.tabix()
        return self.output

    def tabix(self):
        cmd = self.CMD_LINEBREAK.join([
            'tabix',
            '--preset vcf',
            self.output
        ])
        self.call(cmd)


class CopyAndTabixTsvGz(CopyAndTabix):

    seqname_column: int
    start_column: int
    end_column: int

    def main(
            self,
            file: str,
            seqname_column: int,
            start_column: int,
            end_column: int) -> str:

        self.file = file
        self.seqname_column = seqname_column
        self.start_column = start_column
        self.end_column = end_column

        self.copy()
        self.tabix()

        return self.output

    def tabix(self):
        cmd = self.CMD_LINEBREAK.join([
            'tabix',
            f'--sequence {self.seqname_column}',
            f'--begin {self.start_column}',
            f'--end {self.end_column}',
            self.output,
        ])
        self.call(cmd)


class AssertDbnsfpResourceFilenameForVEP(Processor):
    """
    The VEP plugin dbNSFP look for HARD-CODED version string in the filename. Bad software practice!

    After manually checking some possible filenames,
        it seems like it is checking the presence of '4.' or '3.' for versions 4 and 3, respectively
    """

    fname: str

    def main(self, fpath: str):
        self.fname = basename(fpath)
        self.assert_has_version_string()

    def assert_has_version_string(self):
        has_version = ('4.' in self.fname) or ('3.' in self.fname)
        if not has_version:
            error = f'ERROR: dbNSFP resource filename "{self.fname}" does not contain version string "4." or "3."'
            self.logger.info(error)
            raise AssertionError
