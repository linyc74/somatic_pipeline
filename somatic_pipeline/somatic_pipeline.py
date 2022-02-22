from typing import Optional
from .clean_up import CleanUp
from .annotation import SnpEff
from .template import Processor
from .trimming import TrimGalore
from .copy_ref_fa import CopyRefFa
from .parse_vcf import ParseMutect2SnpEffVcf
from .alignment import FactoryIndexAndAlignTumorNormal
from .variant_calling import Mutect2TumorNormalPaired, Mutect2TumorOnly


class SomaticPipeline(Processor):

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    read_aligner: str
    variant_caller: str

    tumor_bam: str
    normal_bam: Optional[str]
    raw_vcf: str
    annotated_vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str],
            read_aligner: str,
            variant_caller: str):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.read_aligner = read_aligner
        self.variant_caller = variant_caller

        self.copy_ref_fa()
        self.trimming()
        self.alignment()
        self.variant_calling()
        self.annotation()
        self.parse_vcf()
        self.clean_up()

    def copy_ref_fa(self):
        self.ref_fa = CopyRefFa(self.settings).main(
            ref_fa=self.ref_fa)

    def trimming(self):
        self.tumor_fq1, self.tumor_fq2 = TrimGalore(self.settings).main(
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2)

        if self.normal_fq1 is not None:
            self.normal_fq1, self.normal_fq2 = TrimGalore(self.settings).main(
                fq1=self.normal_fq1,
                fq2=self.normal_fq2)

    def alignment(self):
        f = FactoryIndexAndAlignTumorNormal().get_callable(
            settings=self.settings,
            read_aligner=self.read_aligner)

        self.tumor_bam, self.normal_bam = f(
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,
            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2)

    def variant_calling(self):
        if self.normal_bam is None:
            self.raw_vcf = Mutect2TumorOnly(self.settings).main(
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam)
        else:
            self.raw_vcf = Mutect2TumorNormalPaired(self.settings).main(
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam)

    def annotation(self):
        self.annotated_vcf = SnpEff(self.settings).main(
            vcf=self.raw_vcf)

    def parse_vcf(self):
        ParseMutect2SnpEffVcf(self.settings).main(
            vcf=self.annotated_vcf)

    def clean_up(self):
        CleanUp(self.settings).main()
