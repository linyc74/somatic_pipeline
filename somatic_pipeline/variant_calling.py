from typing import Optional, Union
from abc import ABC, abstractmethod
from .template import Processor
from .constant import TUMOR, NORMAL


class ProcessInterface(Processor, ABC):

    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]

    tumor_tagged_bam: str
    normal_tagged_bam: Optional[str]
    vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        self.index_ref_fa()

        self.tag_read_groups()
        self.index_bams()

        self.set_vcf()
        self.execute()

        return self.vcf

    def index_ref_fa(self):
        cmd = f'samtools faidx {self.ref_fa}'
        self.call(cmd)

    def tag_read_groups(self):
        self.tumor_tagged_bam = f'{self.workdir}/{TUMOR}-tagged.bam'
        self.__tag_read_group(
            bam_in=self.tumor_bam,
            bam_out=self.tumor_tagged_bam,
            name=TUMOR)

        if self.normal_bam is not None:
            self.normal_tagged_bam = f'{self.workdir}/{NORMAL}-tagged.bam'
            self.__tag_read_group(
                bam_in=self.normal_bam,
                bam_out=self.normal_tagged_bam,
                name=NORMAL)
        else:
            self.normal_tagged_bam = None

    def __tag_read_group(
            self,
            bam_in: str,
            bam_out: str,
            name: str):
        log = f'{self.outdir}/gatk-AddOrReplaceReadGroups-{name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk AddOrReplaceReadGroups',
            f'--INPUT {bam_in}',
            f'--OUTPUT {bam_out}',
            f'--RGLB lib1',
            '--RGPL ILLUMINA',
            '--RGPU unit1',
            f'--RGSM {name}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)

    def index_bams(self):
        self.call(f'samtools index {self.tumor_tagged_bam}')
        if self.normal_tagged_bam is not None:
            self.call(f'samtools index {self.normal_tagged_bam}')

    def set_vcf(self):
        self.vcf = f'{self.workdir}/raw.vcf'

    @abstractmethod
    def execute(self):
        return


class Mutect2Base(ProcessInterface):

    def execute(self):
        return

    def create_sequence_dictionary(self):
        log = f'{self.outdir}/gatk-CreateSequenceDictionary.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk CreateSequenceDictionary',
            f'-R {self.ref_fa}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Mutect2TNPaired(Mutect2Base):

    def execute(self):
        self.create_sequence_dictionary()
        self.mutect2()

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_tagged_bam}',
            f'--input {self.normal_tagged_bam}',
            f'--tumor-sample {TUMOR}',
            f'--normal-sample {NORMAL}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Mutect2TumorOnly(Mutect2Base):

    def execute(self):
        assert self.normal_bam is None
        self.create_sequence_dictionary()
        self.mutect2()

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_tagged_bam}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


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
            self.tumor_tagged_bam,
            self.normal_tagged_bam,
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
            (self.tumor_tagged_bam, self.tumor_pileup, TUMOR),
            (self.normal_tagged_bam, self.normal_pileup, NORMAL),
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
    TUMOR_ONLY_PROCESS_DICT = {
        MUTECT2: Mutect2TumorOnly
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
        Varscan
    ]

    def get_process(
            self,
            variant_caller: str,
            normal_bam: Optional[str]) -> PROCESS_TYPES:

        assert variant_caller in [self.MUTECT2, self.MUSE, self.VARSCAN]

        if normal_bam is None:
            assert variant_caller == self.MUTECT2, \
                f"variant_caller = '{variant_caller}', which does not support tumor-only mode"
            _class = self.TUMOR_ONLY_PROCESS_DICT[variant_caller]

        else:
            _class = self.TN_PAIRED_PROCESS_DICT[variant_caller]

        return _class(self.settings)


class VariantCalling(Processor):

    variant_caller: str
    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]

    vcf: str

    def main(
            self,
            variant_caller: str,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str]) -> str:

        self.variant_caller = variant_caller
        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        process = Factory(self.settings).get_process(
            variant_caller=self.variant_caller,
            normal_bam=self.normal_bam)

        self.vcf = process.main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam)

        return self.vcf
