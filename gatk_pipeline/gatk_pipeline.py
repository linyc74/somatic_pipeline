from .annotation import SnpEff
from .trimming import TrimGalore
from .variant_calling import Mutect2
from .mapping import BwaIndex, BwaMem
from .template import Processor, Settings
from .parse_vcf import ParseMutect2SnpEffVcf
from .constant import TUMOR, NORMAL


class GATKPipeline(Processor):

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: str
    normal_fq2: str

    tumor_bam: str
    normal_bam: str
    raw_vcf: str
    annotated_vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: str,
            normal_fq2: str):

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

    def trimming(self):
        self.tumor_fq1, self.tumor_fq2 = TrimGalore(self.settings).main(
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2)
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
        self.normal_bam = BwaMem(self.settings).main(
            index=index,
            fq1=self.normal_fq1,
            fq2=self.normal_fq2,
            sample_name=NORMAL)

    def variant_calling(self):
        self.raw_vcf = Mutect2(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam)

    def annotation(self):
        self.annotated_vcf = SnpEff(self.settings).main(
            vcf=self.raw_vcf)

    def parse_vcf(self):
        ParseMutect2SnpEffVcf(self.settings).main(
            vcf=self.annotated_vcf)
