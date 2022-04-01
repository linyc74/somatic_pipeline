from typing import Optional, Tuple
from .bqsr import BQSR
from .cnv import ComputeCNV
from .vcf2maf import Vcf2Maf
from .clean_up import CleanUp
from .annotation import SnpEff
from .template import Processor
from .trimming import TrimGalore
from .alignment import Alignment
from .copy_ref_fa import CopyRefFa
from .parse_vcf import ParseSnpEffVcf
from .variant_calling import VariantCalling
from .mark_duplicates import MarkDuplicates


class SomaticPipeline(Processor):

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    read_aligner: str
    variant_caller: str
    exome_target_bed: Optional[str]
    cnvkit_annotate_txt: Optional[str]
    panel_of_normal_vcf: Optional[str]
    bqsr_known_variant_vcf: Optional[str]
    discard_bam: bool
    skip_variant_calling: bool
    skip_cnv: bool

    tumor_bam: str
    normal_bam: Optional[str]

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str],
            read_aligner: str,
            variant_caller: str,
            exome_target_bed: Optional[str],
            cnvkit_annotate_txt: Optional[str],
            panel_of_normal_vcf: Optional[str],
            bqsr_known_variant_vcf: Optional[str],
            discard_bam: bool,
            skip_variant_calling: bool,
            skip_cnv: bool):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.read_aligner = read_aligner
        self.variant_caller = variant_caller
        self.exome_target_bed = exome_target_bed
        self.cnvkit_annotate_txt = cnvkit_annotate_txt
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.bqsr_known_variant_vcf = bqsr_known_variant_vcf
        self.discard_bam = discard_bam
        self.skip_variant_calling = skip_variant_calling
        self.skip_cnv = skip_cnv

        self.copy_ref_fa()
        self.trimming()
        self.alignment_preprocessing_trunk()
        self.variant_calling_trunk()
        self.compute_cnv()
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

    def alignment_preprocessing_trunk(self):
        AlignmentPreprocessingTrunk(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,
            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2,
            read_aligner=self.read_aligner,
            bqsr_known_variant_vcf=self.bqsr_known_variant_vcf,
            discard_bam=self.discard_bam)

    def variant_calling_trunk(self):
        if self.skip_variant_calling:
            self.logger.info(f'Skip all variant calling tasks')
        else:
            VariantCallingTrunk(self.settings).main(
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                variant_caller=self.variant_caller,
                panel_of_normal_vcf=self.panel_of_normal_vcf)

    def compute_cnv(self):
        if self.skip_cnv:
            self.logger.info(f'Skip all CNV tasks')
        else:
            ComputeCNV(self.settings).main(
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                exome_target_bed=self.exome_target_bed,
                annotate_txt=self.cnvkit_annotate_txt)

    def clean_up(self):
        CleanUp(self.settings).main()


class AlignmentPreprocessingTrunk(Processor):

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    read_aligner: str
    bqsr_known_variant_vcf: Optional[str]
    discard_bam: bool

    tumor_bam: str
    normal_bam: Optional[str]

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str],
            read_aligner: str,
            bqsr_known_variant_vcf: Optional[str],
            discard_bam: bool) -> Tuple[str, str]:

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.read_aligner = read_aligner
        self.bqsr_known_variant_vcf = bqsr_known_variant_vcf
        self.discard_bam = discard_bam

        self.alignment()
        self.mark_duplicates()
        self.bqsr()

        return self.tumor_bam, self.normal_bam

    def alignment(self):
        self.tumor_bam, self.normal_bam = Alignment(self.settings).main(
            read_aligner=self.read_aligner,
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,
            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2,
            discard_bam=self.discard_bam)

    def mark_duplicates(self):
        self.tumor_bam, self.normal_bam = MarkDuplicates(self.settings).main(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam)

    def bqsr(self):
        self.tumor_bam, self.normal_bam = BQSR(self.settings).main(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            ref_fa=self.ref_fa,
            known_variant_vcf=self.bqsr_known_variant_vcf)


class VariantCallingTrunk(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]
    variant_caller: str
    panel_of_normal_vcf: Optional[str]

    raw_vcf: str
    annotated_vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            variant_caller: str,
            panel_of_normal_vcf: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_caller = variant_caller
        self.panel_of_normal_vcf = panel_of_normal_vcf

        self.variant_calling()
        self.annotation()
        self.parse_vcf()
        self.vcf_2_maf()

    def variant_calling(self):
        self.raw_vcf = VariantCalling(self.settings).main(
            variant_caller=self.variant_caller,
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=self.panel_of_normal_vcf)

    def annotation(self):
        self.annotated_vcf = SnpEff(self.settings).main(
            vcf=self.raw_vcf)

    def parse_vcf(self):
        ParseSnpEffVcf(self.settings).main(
            vcf=self.annotated_vcf)

    def vcf_2_maf(self):
        Vcf2Maf(self.settings).main(
            annotated_vcf=self.annotated_vcf,
            ref_fa=self.ref_fa,
            variant_caller=self.variant_caller)
