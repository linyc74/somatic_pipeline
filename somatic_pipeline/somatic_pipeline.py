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
from .variant_picking import VariantPicking


class SomaticPipeline(Processor):

    # required
    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str

    # optional
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]

    # preprocessing
    clip_r1_5_prime: int
    clip_r2_5_prime: int
    read_aligner: str
    skip_mark_duplicates: bool
    bqsr_known_variant_vcf: Optional[str]
    discard_bam: bool
    tumor_bam: str
    normal_bam: Optional[str]

    # variant calling
    variant_callers: List[str]
    skip_variant_calling: bool
    call_region_bed: Optional[str]
    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    # variant filtering
    variant_flagging_criteria: Optional[str]
    variant_removal_flags: List[str]

    # variant picking
    min_snv_callers: int
    min_indel_callers: int

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

            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            read_aligner: str,
            skip_mark_duplicates: bool,
            bqsr_known_variant_vcf: Optional[str],
            discard_bam: bool,

            variant_callers: List[str],
            skip_variant_calling: bool,
            call_region_bed: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],

            variant_flagging_criteria: Optional[str],
            variant_removal_flags: List[str],

            min_snv_callers: int,
            min_indel_callers: int,

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

        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime
        self.read_aligner = read_aligner
        self.skip_mark_duplicates = skip_mark_duplicates
        self.bqsr_known_variant_vcf = bqsr_known_variant_vcf
        self.discard_bam = discard_bam

        self.variant_callers = variant_callers
        self.skip_variant_calling = skip_variant_calling
        self.call_region_bed = call_region_bed
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.variant_flagging_criteria = variant_flagging_criteria
        self.variant_removal_flags = variant_removal_flags

        self.min_snv_callers = min_snv_callers
        self.min_indel_callers = min_indel_callers

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
            fq2=self.tumor_fq2,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)

        if self.normal_fq1 is not None:
            self.normal_fq1, self.normal_fq2 = TrimGalore(self.settings).main(
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                clip_r1_5_prime=self.clip_r1_5_prime,
                clip_r2_5_prime=self.clip_r2_5_prime)

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
                variant_callers=self.variant_callers,
                call_region_bed=self.call_region_bed,
                panel_of_normal_vcf=self.panel_of_normal_vcf,
                germline_resource_vcf=self.germline_resource_vcf,
                variant_flagging_criteria=self.variant_flagging_criteria,
                variant_removal_flags=self.variant_removal_flags,
                min_snv_callers=self.min_snv_callers,
                min_indel_callers=self.min_indel_callers,
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

    variant_callers: List[str]
    call_region_bed: Optional[str]
    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    variant_flagging_criteria: Optional[str]
    variant_removal_flags: List[str]

    min_snv_callers: int
    min_indel_callers: int

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

    vcfs: List[str]
    vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],

            variant_callers: List[str],
            call_region_bed: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],

            variant_flagging_criteria: Optional[str],
            variant_removal_flags: List[str],

            min_snv_callers: int,
            min_indel_callers: int,

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

        self.variant_callers = variant_callers
        self.call_region_bed = call_region_bed
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.variant_flagging_criteria = variant_flagging_criteria
        self.variant_removal_flags = variant_removal_flags

        self.min_snv_callers = min_snv_callers
        self.min_indel_callers = min_indel_callers

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
        self.variant_picking()
        self.annotation()
        self.parse_vcf()
        self.vcf_2_maf()
        self.compress_index_vcf()

    def variant_calling(self):
        self.vcfs = VariantCalling(self.settings).main(
            variant_callers=self.variant_callers,
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            call_region_bed=self.call_region_bed,
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            germline_resource_vcf=self.germline_resource_vcf,
            variant_flagging_criteria=self.variant_flagging_criteria,
            variant_removal_flags=self.variant_removal_flags)

    def variant_picking(self):
        self.vcf = VariantPicking(self.settings).main(
            ref_fa=self.ref_fa,
            vcfs=self.vcfs,
            min_snv_callers=self.min_snv_callers,
            min_indel_callers=self.min_indel_callers)

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

    def parse_vcf(self):
        for vcf in self.vcfs + [self.vcf]:
            ParseVcf(self.settings).main(
                vcf=vcf,
                dstdir=None  # same dir of input vcf
            )

    def vcf_2_maf(self):
        for vcf in self.vcfs + [self.vcf]:
            Vcf2Maf(self.settings).main(
                vcf=vcf,
                ref_fa=self.ref_fa,
                dstdir=None  # same dir of input vcf
            )

    def compress_index_vcf(self):
        for vcf in self.vcfs + [self.vcf]:
            BgzipIndex(self.settings).main(vcf=vcf, keep=False)
