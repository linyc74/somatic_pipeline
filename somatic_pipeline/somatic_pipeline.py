from typing import Optional, Tuple, List
from .bqsr import BQSR
from .cnv import ComputeCNV
from .vcf2maf import Vcf2Maf
from .clean_up import CleanUp
from .parse_vcf import ParseVcf
from .template import Processor
from .trimming import TrimGalore
from .alignment import Alignment
from .annotation import Annotation
from .copy_ref_fa import CopyRefFa
from .map_stats import MappingStats
from .index_files import BgzipIndex
from .variant_calling import VariantCalling
from .mark_duplicates import MarkDuplicates
from .variant_filtering import Mutect2VariantFiltering


class SomaticPipeline(Processor):

    # required
    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str

    # optional
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]

    # preprocessing
    read_aligner: str
    skip_mark_duplicates: bool
    bqsr_known_variant_vcf: Optional[str]
    discard_bam: bool
    tumor_bam: str
    normal_bam: Optional[str]

    # variant calling
    variant_caller: str
    skip_variant_calling: bool
    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    # mutect2 variant filtering
    filter_mutect2_variants: bool
    variant_removal_flags: List[str]

    # annotation
    annotator: str
    skip_variant_annotation: bool
    vep_db_tar_gz: Optional[str]
    vep_db_type: str
    vep_buffer_size: int
    dbnsfp_resource: Optional[str]
    cadd_resource: Optional[str]
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]

    # cnv
    skip_cnv: bool
    cnvkit_annotate_txt: Optional[str]
    exome_target_bed: Optional[str]

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,

            normal_fq1: Optional[str],
            normal_fq2: Optional[str],

            read_aligner: str,
            skip_mark_duplicates: bool,
            bqsr_known_variant_vcf: Optional[str],
            discard_bam: bool,

            variant_caller: str,
            skip_variant_calling: bool,
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],

            filter_mutect2_variants: bool,
            variant_removal_flags: List[str],

            annotator: str,
            skip_variant_annotation: bool,
            vep_db_tar_gz: Optional[str],
            vep_db_type: str,
            vep_buffer_size: int,
            dbnsfp_resource: Optional[str],
            cadd_resource: Optional[str],
            clinvar_vcf_gz: Optional[str],
            dbsnp_vcf_gz: Optional[str],
            snpsift_dbnsfp_txt_gz: Optional[str],

            skip_cnv: bool,
            exome_target_bed: Optional[str],
            cnvkit_annotate_txt: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2

        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2

        self.read_aligner = read_aligner
        self.skip_mark_duplicates = skip_mark_duplicates
        self.bqsr_known_variant_vcf = bqsr_known_variant_vcf
        self.discard_bam = discard_bam

        self.variant_caller = variant_caller
        self.skip_variant_calling = skip_variant_calling
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.filter_mutect2_variants = filter_mutect2_variants
        self.variant_removal_flags = variant_removal_flags

        self.annotator = annotator
        self.skip_variant_annotation = skip_variant_annotation
        self.vep_db_tar_gz = vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.vep_buffer_size = vep_buffer_size
        self.dbnsfp_resource = dbnsfp_resource
        self.cadd_resource = cadd_resource
        self.clinvar_vcf_gz = clinvar_vcf_gz
        self.dbsnp_vcf_gz = dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = snpsift_dbnsfp_txt_gz

        self.skip_cnv = skip_cnv
        self.exome_target_bed = exome_target_bed
        self.cnvkit_annotate_txt = cnvkit_annotate_txt

        self.copy_ref_fa()
        self.trimming()
        self.alignment_preprocessing_workflow()
        self.variant_calling_workflow()
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

    def alignment_preprocessing_workflow(self):
        self.tumor_bam, self.normal_bam = AlignmentPreprocessingWorkflow(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,
            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2,
            read_aligner=self.read_aligner,
            bqsr_known_variant_vcf=self.bqsr_known_variant_vcf,
            discard_bam=self.discard_bam,
            skip_mark_duplicates=self.skip_mark_duplicates)

    def variant_calling_workflow(self):
        if self.skip_variant_calling:
            self.logger.info(f'Skip all variant calling tasks')
        else:
            VariantCallingWorkflow(self.settings).main(
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                variant_caller=self.variant_caller,
                panel_of_normal_vcf=self.panel_of_normal_vcf,
                germline_resource_vcf=self.germline_resource_vcf,
                filter_mutect2_variants=self.filter_mutect2_variants,
                variant_removal_flags=self.variant_removal_flags,
                annotator=self.annotator,
                skip_variant_annotation=self.skip_variant_annotation,
                clinvar_vcf_gz=self.clinvar_vcf_gz,
                dbsnp_vcf_gz=self.dbsnp_vcf_gz,
                snpsift_dbnsfp_txt_gz=self.snpsift_dbnsfp_txt_gz,
                vep_db_tar_gz=self.vep_db_tar_gz,
                vep_db_type=self.vep_db_type,
                vep_buffer_size=self.vep_buffer_size,
                cadd_resource=self.cadd_resource,
                dbnsfp_resource=self.dbnsfp_resource)

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


class AlignmentPreprocessingWorkflow(Processor):

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    read_aligner: str
    bqsr_known_variant_vcf: Optional[str]
    discard_bam: bool
    skip_mark_duplicates: bool

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
            discard_bam: bool,
            skip_mark_duplicates: bool) -> Tuple[str, str]:

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.read_aligner = read_aligner
        self.bqsr_known_variant_vcf = bqsr_known_variant_vcf
        self.discard_bam = discard_bam
        self.skip_mark_duplicates = skip_mark_duplicates

        self.alignment()
        self.mark_duplicates()
        self.bqsr()
        self.mapping_stats()

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
        if self.skip_mark_duplicates:
            self.logger.info(f'Skip mark PCR duplicates')
            return
        self.tumor_bam, self.normal_bam = MarkDuplicates(self.settings).main(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam)

    def bqsr(self):
        if self.bqsr_known_variant_vcf is None:
            self.logger.info(f'BQSR known variant VCF not provided, skip BQSR')
            return
        self.tumor_bam, self.normal_bam = BQSR(self.settings).main(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            ref_fa=self.ref_fa,
            known_variant_vcf=self.bqsr_known_variant_vcf)

    def mapping_stats(self):
        MappingStats(self.settings).main(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam)


class VariantCallingWorkflow(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]

    variant_caller: str
    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    filter_mutect2_variants: bool
    variant_removal_flags: List[str]

    annotator: str
    skip_variant_annotation: bool
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]
    vep_db_tar_gz: Optional[str]
    vep_db_type: str
    vep_buffer_size: int
    cadd_resource: Optional[str]
    dbnsfp_resource: Optional[str]

    vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],

            variant_caller: str,
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],

            filter_mutect2_variants: bool,
            variant_removal_flags: List[str],

            annotator: str,
            skip_variant_annotation: bool,
            clinvar_vcf_gz: Optional[str],
            dbsnp_vcf_gz: Optional[str],
            snpsift_dbnsfp_txt_gz: Optional[str],
            vep_db_tar_gz: Optional[str],
            vep_db_type: str,
            vep_buffer_size: int,
            cadd_resource: Optional[str],
            dbnsfp_resource: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        self.variant_caller = variant_caller
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.filter_mutect2_variants = filter_mutect2_variants
        self.variant_removal_flags = variant_removal_flags

        self.annotator = annotator
        self.skip_variant_annotation = skip_variant_annotation
        self.clinvar_vcf_gz = clinvar_vcf_gz
        self.dbsnp_vcf_gz = dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = snpsift_dbnsfp_txt_gz
        self.vep_db_tar_gz = vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.vep_buffer_size = vep_buffer_size
        self.cadd_resource = cadd_resource
        self.dbnsfp_resource = dbnsfp_resource

        self.variant_calling()
        self.mutect2_variant_filtering()
        self.annotation()
        self.move_vcf_to_outdir()
        self.parse_vcf()
        self.vcf_2_maf()
        self.compress_index_vcf()

    def variant_calling(self):
        self.vcf = VariantCalling(self.settings).main(
            variant_caller=self.variant_caller,
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            germline_resource_vcf=self.germline_resource_vcf)

    def mutect2_variant_filtering(self):
        if self.filter_mutect2_variants:
            self.vcf = Mutect2VariantFiltering(self.settings).main(
                vcf=self.vcf,
                ref_fa=self.ref_fa,
                variant_caller=self.variant_caller,
                variant_removal_flags=self.variant_removal_flags)

    def annotation(self):
        if not self.skip_variant_annotation:
            self.vcf = Annotation(self.settings).main(
                annotator=self.annotator,
                vcf=self.vcf,
                ref_fa=self.ref_fa,
                clinvar_vcf_gz=self.clinvar_vcf_gz,
                dbsnp_vcf_gz=self.dbsnp_vcf_gz,
                snpsift_dbnsfp_txt_gz=self.snpsift_dbnsfp_txt_gz,
                vep_db_tar_gz=self.vep_db_tar_gz,
                vep_db_type=self.vep_db_type,
                vep_buffer_size=self.vep_buffer_size,
                cadd_resource=self.cadd_resource,
                dbnsfp_resource=self.dbnsfp_resource)

    def move_vcf_to_outdir(self):
        dst = f'{self.outdir}/variants.vcf'
        self.call(f'mv {self.vcf} {dst}')
        self.vcf = dst

    def parse_vcf(self):
        ParseVcf(self.settings).main(vcf=self.vcf)

    def vcf_2_maf(self):
        Vcf2Maf(self.settings).main(
            annotated_vcf=self.vcf,
            ref_fa=self.ref_fa,
            variant_caller=self.variant_caller)

    def compress_index_vcf(self):
        BgzipIndex(self.settings).main(vcf=self.vcf, keep=False)
