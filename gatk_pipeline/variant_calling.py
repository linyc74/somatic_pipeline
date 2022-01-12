from os.path import basename
from .constant import CMD_LINEBREAK
from .constant import TUMOR, NORMAL
from .template import Processor, Settings


class Mutect2(Processor):

    ref_fa: str
    copied_ref_fa: str
    vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def prepare_ref_fa(self):
        self.copy_ref_fa()
        self.index_ref_fa()

    def copy_ref_fa(self):
        self.copied_ref_fa = f'{self.workdir}/{basename(self.ref_fa)}'
        cmd = f'cp {self.ref_fa} {self.copied_ref_fa}'
        self.call(cmd)

    def index_ref_fa(self):
        cmd = f'samtools faidx {self.copied_ref_fa}'
        self.call(cmd)

        log = f'{self.outdir}/gatk_CreateSequenceDictionary.log'
        cmd = CMD_LINEBREAK.join([
            'gatk CreateSequenceDictionary',
            f'-R {self.copied_ref_fa}',
            f'1> {log} 2> {log}'
        ])
        self.call(cmd)

    def tag_read_group(
            self,
            bam_in: str,
            bam_out: str,
            name: str):
        log = f'{self.outdir}/gatk_AddOrReplaceReadGroups_{name}.log'
        cmd = CMD_LINEBREAK.join([
            'gatk AddOrReplaceReadGroups',
            f'--INPUT {bam_in}',
            f'--OUTPUT {bam_out}',
            f'--RGLB lib1',
            '--RGPL ILLUMINA',
            '--RGPU unit1',
            f'--RGSM {name}',
            f'1> {log} 2> {log}',
        ])
        self.call(cmd)

    def index_bam(self, bam: str):
        self.call(f'samtools index {bam}')


class Mutect2TumorNormalPaired(Mutect2):

    tumor_bam: str
    normal_bam: str

    tumor_tagged_bam: str
    normal_tagged_bam: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        self.prepare_ref_fa()
        self.tag_tumor_and_normal_read_groups()
        self.index_tumor_normal_bams()
        self.mutect2()

        return self.vcf

    def tag_tumor_and_normal_read_groups(self):
        self.normal_tagged_bam = f'{self.workdir}/{NORMAL}_read_group_tagged.bam'
        self.tumor_tagged_bam = f'{self.workdir}/{TUMOR}_read_group_tagged.bam'

        for bam_in, bam_out, name in [
            (self.normal_bam, self.normal_tagged_bam, NORMAL),
            (self.tumor_bam, self.tumor_tagged_bam, TUMOR),
        ]:
            self.tag_read_group(bam_in=bam_in, bam_out=bam_out, name=name)

    def index_tumor_normal_bams(self):
        for bam in [self.normal_tagged_bam, self.tumor_tagged_bam]:
            self.index_bam(bam=bam)

    def mutect2(self):
        self.vcf = f'{self.workdir}/raw.vcf'
        log = f'{self.outdir}/gatk_Mutect2.log'
        cmd = CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.copied_ref_fa}',
            f'--input {self.tumor_tagged_bam}',
            f'--input {self.normal_tagged_bam}',
            f'--normal-sample {NORMAL}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'1> {log} 2> {log}',
        ])
        self.call(cmd)


class Mutect2TumorOnly(Mutect2):

    tumor_bam: str
    tumor_tagged_bam: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam

        self.prepare_ref_fa()
        self.tag_tumor_read_group()
        self.index_tumor_bam()
        self.mutect2()

        return self.vcf

    def tag_tumor_read_group(self):
        self.tumor_tagged_bam = f'{self.workdir}/{TUMOR}_read_group_tagged.bam'
        self.tag_read_group(
            bam_in=self.tumor_bam,
            bam_out=self.tumor_tagged_bam,
            name=TUMOR)

    def index_tumor_bam(self):
        self.index_bam(bam=self.tumor_tagged_bam)

    def mutect2(self):
        self.vcf = f'{self.workdir}/raw.vcf'
        log = f'{self.outdir}/gatk_Mutect2.log'
        cmd = CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.copied_ref_fa}',
            f'--input {self.tumor_tagged_bam}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'1> {log} 2> {log}',
        ])
        self.call(cmd)
