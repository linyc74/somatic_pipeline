from os.path import basename
from .constant import CMD_LINEBREAK
from .template import Processor, Settings


class Mutect2(Processor):

    TUMOR_NAME = 'tumor'
    NORMAL_NAME = 'normal'

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    vcf: str

    copied_ref_fa: str
    tumor_tagged_bam: str
    normal_tagged_bam: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        self.prepare_ref_fa()
        self.tag_read_group()
        self.index_bam()
        self.mutect2()

        return self.vcf

    def prepare_ref_fa(self):
        self.__copy_ref_fa()
        self.__index_ref_fa()

    def __copy_ref_fa(self):
        self.copied_ref_fa = f'{self.workdir}/{basename(self.ref_fa)}'
        cmd = f'cp {self.ref_fa} {self.copied_ref_fa}'
        self.call(cmd)

    def __index_ref_fa(self):
        cmd = f'samtools faidx {self.copied_ref_fa}'
        self.call(cmd)
        cmd = f'gatk CreateSequenceDictionary -R {self.copied_ref_fa}'
        self.call(cmd)

    def tag_read_group(self):
        self.normal_tagged_bam = f'{self.workdir}/{self.NORMAL_NAME}_read_group_tagged.bam'
        self.tumor_tagged_bam = f'{self.workdir}/{self.TUMOR_NAME}_read_group_tagged.bam'

        for bam_in, bam_out, name in [
            (self.normal_bam, self.normal_tagged_bam, self.NORMAL_NAME),
            (self.tumor_bam, self.tumor_tagged_bam, self.TUMOR_NAME),
        ]:
            cmd = CMD_LINEBREAK.join([
                'gatk AddOrReplaceReadGroups',
                f'--INPUT {bam_in}',
                f'--OUTPUT {bam_out}',
                f'--RGLB lib1',
                '--RGPL ILLUMINA',
                '--RGPU unit1',
                f'--RGSM {name}'
            ])
            self.call(cmd)

    def index_bam(self):
        for bam in [self.normal_tagged_bam, self.tumor_tagged_bam]:
            cmd = f'samtools index {bam}'
            self.call(cmd)

    def mutect2(self):
        self.vcf = f'{self.outdir}/raw.vcf'
        cmd = CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.copied_ref_fa}',
            f'--input {self.tumor_tagged_bam}',
            f'--input {self.normal_tagged_bam}',
            f'--normal-sample {self.NORMAL_NAME}',
            f'--output {self.vcf}',
        ])
        self.call(cmd)
