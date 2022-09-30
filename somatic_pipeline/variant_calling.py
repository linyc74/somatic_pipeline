from abc import ABC, abstractmethod
from typing import Optional, Union, List
from .tools import edit_fpath
from .template import Processor
from .constant import TUMOR, NORMAL
from .index_files import SamtoolsIndexFa, GATKIndexVcf, GATKCreateSequenceDictionary, SamtoolsIndexBam


class ProcessInterface(Processor, ABC):

    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]
    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.index_ref_fa()
        self.index_bams()
        self.set_vcf()
        self.execute()

        return self.vcf

    def index_ref_fa(self):
        SamtoolsIndexFa(self.settings).main(fa=self.ref_fa)

    def index_bams(self):
        SamtoolsIndexBam(self.settings).main(bam=self.tumor_bam)
        if self.normal_bam is not None:
            SamtoolsIndexBam(self.settings).main(bam=self.normal_bam)

    def set_vcf(self):
        self.vcf = f'{self.workdir}/raw.vcf'

    @abstractmethod
    def execute(self):
        return


class GATKProcessInterface(ProcessInterface):

    pon_args: List[str]
    germline_resource_args: List[str]

    def execute(self):
        return

    def create_sequence_dictionary(self):
        GATKCreateSequenceDictionary(self.settings).main(ref_fa=self.ref_fa)

    def prepare_resource_vcfs(self):
        self.pon_args, self.germline_resource_args = PrepareResourceVcfs(self.settings).main(
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            germline_resource_vcf=self.germline_resource_vcf)


class PrepareResourceVcfs(Processor):

    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    pon_args: List[str]
    germline_resource_args: List[str]

    def main(
            self,
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str]):

        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.prepare_panel_of_normal()
        self.prepare_germline_resource()

        return self.pon_args, self.germline_resource_args

    def prepare_panel_of_normal(self):
        if self.panel_of_normal_vcf is None:
            self.pon_args = []
        else:
            self.panel_of_normal_vcf = self.__copy(self.panel_of_normal_vcf)
            self.__index(self.panel_of_normal_vcf)
            self.pon_args = [f'--panel-of-normals {self.panel_of_normal_vcf}']

    def prepare_germline_resource(self):
        if self.germline_resource_vcf is None:
            self.germline_resource_args = []
        else:
            self.germline_resource_vcf = self.__copy(self.germline_resource_vcf)
            self.__index(self.germline_resource_vcf)
            self.germline_resource_args = [f'--germline-resource {self.germline_resource_vcf}']

    def __copy(self, vcf: str) -> str:
        dst = edit_fpath(fpath=vcf, dstdir=self.workdir)
        self.call(f'cp {vcf} {dst}')
        return dst

    def __index(self, vcf: str):
        GATKIndexVcf(self.settings).main(vcf=vcf)


class Mutect2TNPaired(GATKProcessInterface):

    def execute(self):
        self.create_sequence_dictionary()
        self.prepare_resource_vcfs()
        self.mutect2()

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_bam}',
            f'--input {self.normal_bam}',
            f'--tumor-sample {TUMOR}',
            f'--normal-sample {NORMAL}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
        ] + self.pon_args + self.germline_resource_args + [
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Mutect2TumorOnly(GATKProcessInterface):

    def execute(self):
        assert self.normal_bam is None
        self.create_sequence_dictionary()
        self.prepare_resource_vcfs()
        self.mutect2()

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_bam}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
        ] + self.pon_args + self.germline_resource_args + [
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class HaplotypeCaller(GATKProcessInterface):

    def execute(self):
        assert self.normal_bam is None
        self.create_sequence_dictionary()
        self.haplotype_caller()

    def haplotype_caller(self):
        log = f'{self.outdir}/gatk-HaplotypeCaller.log'
        args = [
            'gatk HaplotypeCaller',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_bam}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))


class Muse(ProcessInterface):

    WES_OR_WGS = '-E'  # '-E' for WES; '-G' for WGS

    call_result_txt: str

    def execute(self):
        self.muse_call()
        self.muse_sump()

    def muse_call(self):
        log = f'{self.outdir}/MuSE-call.log'
        output = f'{self.workdir}/MuSE-call-result'
        cmd = self.CMD_LINEBREAK.join([
            'MuSE call',
            f'-f {self.ref_fa}',
            f'-O {output}',
            self.tumor_bam,
            self.normal_bam,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)
        self.call_result_txt = output + '.MuSE.txt'

    def muse_sump(self):
        log = f'{self.outdir}/MuSE-sump.log'
        cmd = self.CMD_LINEBREAK.join([
            'MuSE sump',
            f'-I {self.call_result_txt}',
            self.WES_OR_WGS,
            f'-O {self.vcf}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Varscan(ProcessInterface):

    normal_pileup: str
    tumor_pileup: str
    snp_vcf: str
    indel_vcf: str

    def execute(self):
        self.samtools_mpileup()
        self.varscan_somatic()
        self.compress_vcfs()
        self.index_vcfs()
        self.concat_snp_indel_vcfs()

    def samtools_mpileup(self):
        self.tumor_pileup = f'{self.workdir}/{TUMOR}-pileup'
        self.normal_pileup = f'{self.workdir}/{NORMAL}-pileup'
        for bam, output, name in [
            (self.tumor_bam, self.tumor_pileup, TUMOR),
            (self.normal_bam, self.normal_pileup, NORMAL),
        ]:
            self.__mpileup(bam=bam, output=output, name=name)

    def __mpileup(self, bam: str, output: str, name: str):
        log = f'{self.outdir}/samtools-mpileup-{name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'samtools mpileup',
            f'--fasta-ref {self.ref_fa}',
            f'--output {output}',
            bam,
            f'2> {log}',
        ])
        self.call(cmd)

    def varscan_somatic(self):
        self.snp_vcf = f'{self.workdir}/varscan-snp.vcf'
        self.indel_vcf = f'{self.workdir}/varscan-indel.vcf'
        log = f'{self.outdir}/varscan-somatic.log'
        cmd = self.CMD_LINEBREAK.join([
            'varscan somatic',
            self.normal_pileup,
            self.tumor_pileup,
            f'--output-snp {self.snp_vcf}',
            f'--output-indel {self.indel_vcf}',
            '--strand-filter 1',  # removes variants with >90% strand bias
            '--output-vcf 1',  # output VCF instead of VarScan native format
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)

    def compress_vcfs(self):
        for vcf in [self.snp_vcf, self.indel_vcf]:
            self.call(f'bgzip --force --threads {self.threads} {vcf}')
        self.snp_vcf += '.gz'
        self.indel_vcf += '.gz'

    def index_vcfs(self):
        for vcf in [self.snp_vcf, self.indel_vcf]:
            self.call(f'bcftools index --force --threads {self.threads} {vcf}')

    def concat_snp_indel_vcfs(self):
        log = f'{self.outdir}/bcftools-merge.log'
        cmd = self.CMD_LINEBREAK.join([
            'bcftools concat',
            '--allow-overlaps',  # first coordinate of the next file can precede last record of the current file
            f'--threads {self.threads}',
            f'--output-type v',  # uncompressed VCF
            f'--output {self.vcf}',
            self.snp_vcf,
            self.indel_vcf,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Factory(Processor):

    MUTECT2 = 'mutect2'
    MUSE = 'muse'
    VARSCAN = 'varscan'
    HAPLOTYPE_CALLER = 'haplotype-caller'
    TUMOR_ONLY_PROCESS_DICT = {
        MUTECT2: Mutect2TumorOnly,
        HAPLOTYPE_CALLER: HaplotypeCaller
    }
    TN_PAIRED_PROCESS_DICT = {
        MUTECT2: Mutect2TNPaired,
        MUSE: Muse,
        VARSCAN: Varscan
    }
    PROCESS_TYPES = Union[
        Mutect2TNPaired,
        Mutect2TumorOnly,
        Muse,
        Varscan,
        HaplotypeCaller,
    ]

    def get_process(
            self,
            variant_caller: str,
            normal_bam: Optional[str]) -> PROCESS_TYPES:

        assert variant_caller in [self.MUTECT2, self.MUSE, self.VARSCAN, self.HAPLOTYPE_CALLER]

        if normal_bam is None:
            assert variant_caller in [self.MUTECT2, self.HAPLOTYPE_CALLER], \
                f"variant_caller = '{variant_caller}', which does not support tumor-only mode"
            _class = self.TUMOR_ONLY_PROCESS_DICT[variant_caller]

        else:
            _class = self.TN_PAIRED_PROCESS_DICT[variant_caller]

        return _class(self.settings)


class VariantCalling(Processor):

    def main(
            self,
            variant_caller: str,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str]) -> str:

        process = Factory(self.settings).get_process(
            variant_caller=variant_caller,
            normal_bam=normal_bam)

        vcf = process.main(
            ref_fa=ref_fa,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            panel_of_normal_vcf=panel_of_normal_vcf,
            germline_resource_vcf=germline_resource_vcf)

        return vcf
