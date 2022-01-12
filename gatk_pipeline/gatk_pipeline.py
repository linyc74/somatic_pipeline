from typing import Optional
from .clean_up import CleanUp
from .annotation import SnpEff
from .trimming import TrimGalore
from .constant import TUMOR, NORMAL
from .mapping import BwaIndex, BwaMem
from .template import Processor, Settings
from .parse_vcf import ParseMutect2SnpEffVcf
from .variant_calling import Mutect2TumorNormalPaired, Mutect2TumorOnly


class GATKPipeline(Processor):

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]

    tumor_bam: str
    normal_bam: Optional[str]
    raw_vcf: str
    annotated_vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2

        self.trimming()
        self.mapping()
        self.variant_calling()
        self.annotation()
        self.parse_vcf()
        self.clean_up()

    def trimming(self):
        self.tumor_fq1, self.tumor_fq2 = TrimGalore(self.settings).main(
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2)

        if self.normal_fq1 is None:
            return
        self.normal_fq1, self.normal_fq2 = TrimGalore(self.settings).main(
            fq1=self.normal_fq1,
            fq2=self.normal_fq2)

    def mapping(self):
        index = BwaIndex(self.settings).main(
            fna=self.ref_fa)

        self.tumor_bam = BwaMem(self.settings).main(
            index=index,
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2,
            sample_name=TUMOR)

        if self.normal_fq1 is None:
            self.normal_bam = None
        else:
            self.normal_bam = BwaMem(self.settings).main(
                index=index,
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                sample_name=NORMAL)

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
