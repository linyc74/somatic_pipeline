import os
from typing import Optional
from .template import Settings
from .tools import get_temp_path
from .somatic_pipeline import SomaticPipeline


class Main:

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str

    normal_fq1: Optional[str]
    normal_fq2: Optional[str]

    read_aligner: str
    skip_mark_duplicates: bool
    bqsr_known_variant_vcf: Optional[str]
    discard_bam: bool

    variant_caller: str
    skip_variant_calling: bool
    panel_of_normal_vcf: Optional[str]

    annotator: str
    vep_db_tar_gz: Optional[str]
    vep_db_type: str
    vep_buffer_size: int
    dbnsfp_resource: Optional[str]
    cadd_resource: Optional[str]
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]

    skip_cnv: bool
    exome_target_bed: Optional[str]
    cnvkit_annotate_txt: Optional[str]

    settings: Settings

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,

            normal_fq1: str,
            normal_fq2: str,
            outdir: str,
            threads: str,
            debug: bool,

            read_aligner: str,
            skip_mark_duplicates: bool,
            bqsr_known_variant_vcf: str,
            discard_bam: bool,

            variant_caller: str,
            skip_variant_calling: bool,
            panel_of_normal_vcf: str,

            annotator: str,
            vep_db_tar_gz: str,
            vep_db_type: str,
            vep_buffer_size: int,
            dbnsfp_resource: str,
            cadd_resource: str,
            clinvar_vcf_gz: str,
            dbsnp_vcf_gz: str,
            snpsift_dbnsfp_txt_gz: str,

            skip_cnv: bool,
            exome_target_bed: str,
            cnvkit_annotate_txt: str):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2

        self.normal_fq1 = None if normal_fq1.lower() == 'none' else normal_fq1
        self.normal_fq2 = None if normal_fq2.lower() == 'none' else normal_fq2

        self.read_aligner = read_aligner
        self.skip_mark_duplicates = skip_mark_duplicates
        self.bqsr_known_variant_vcf = None if bqsr_known_variant_vcf.lower() == 'none' else bqsr_known_variant_vcf
        self.discard_bam = discard_bam

        self.variant_caller = variant_caller
        self.skip_variant_calling = skip_variant_calling
        self.panel_of_normal_vcf = None if panel_of_normal_vcf.lower() == 'none' else panel_of_normal_vcf

        self.annotator = annotator
        self.vep_db_tar_gz = None if vep_db_tar_gz.lower() == 'none' else vep_db_tar_gz
        self.vep_db_type = vep_db_type
        self.vep_buffer_size = vep_buffer_size
        self.dbnsfp_resource = None if dbnsfp_resource.lower() == 'none' else dbnsfp_resource
        self.cadd_resource = None if cadd_resource.lower() == 'none' else cadd_resource
        self.clinvar_vcf_gz = None if clinvar_vcf_gz.lower() == 'none' else clinvar_vcf_gz
        self.dbsnp_vcf_gz = None if dbsnp_vcf_gz.lower() == 'none' else dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = None if snpsift_dbnsfp_txt_gz.lower() == 'none' else snpsift_dbnsfp_txt_gz

        self.skip_cnv = skip_cnv
        self.exome_target_bed = None if exome_target_bed.lower() == 'none' else exome_target_bed
        self.cnvkit_annotate_txt = None if cnvkit_annotate_txt.lower() == 'none' else cnvkit_annotate_txt

        self.settings = Settings(
            workdir=get_temp_path(prefix='./somatic_pipeline_workdir_'),
            outdir=outdir,
            threads=int(threads),
            debug=debug,
            mock=False)

        for d in [self.settings.workdir, self.settings.outdir]:
            os.makedirs(d, exist_ok=True)

        SomaticPipeline(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,

            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2,

            read_aligner=self.read_aligner,
            skip_mark_duplicates=self.skip_mark_duplicates,
            bqsr_known_variant_vcf=self.bqsr_known_variant_vcf,
            discard_bam=self.discard_bam,

            variant_caller=self.variant_caller,
            skip_variant_calling=self.skip_variant_calling,
            panel_of_normal_vcf=self.panel_of_normal_vcf,

            annotator=self.annotator,
            vep_db_tar_gz=self.vep_db_tar_gz,
            vep_db_type=self.vep_db_type,
            vep_buffer_size=self.vep_buffer_size,
            dbnsfp_resource=self.dbnsfp_resource,
            cadd_resource=self.cadd_resource,
            clinvar_vcf_gz=self.clinvar_vcf_gz,
            dbsnp_vcf_gz=self.dbsnp_vcf_gz,
            snpsift_dbnsfp_txt_gz=self.snpsift_dbnsfp_txt_gz,

            skip_cnv=self.skip_cnv,
            exome_target_bed=self.exome_target_bed,
            cnvkit_annotate_txt=self.cnvkit_annotate_txt)
