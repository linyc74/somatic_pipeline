import os
from .template import Settings
from .tools import get_temp_path
from .annotation import Annotation
from .somatic_pipeline import SomaticPipeline


class Run:

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

            clip_r1_5_prime: int,
            clip_r2_5_prime: int,
            read_aligner: str,
            skip_mark_duplicates: bool,
            bqsr_known_variant_vcf: str,
            discard_bam: bool,

            variant_callers: str,
            skip_variant_calling: bool,
            panel_of_normal_vcf: str,
            germline_resource_vcf: str,
            call_region_bed: str,

            variant_flagging_criteria: str,
            variant_removal_flags: str,

            min_snv_callers: int,
            min_indel_callers: int,

            annotator: str,
            skip_variant_annotation: bool,
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

        self.config_settings(outdir=outdir, threads=int(threads), debug=debug)

        SomaticPipeline(self.settings).main(
            ref_fa=ref_fa,
            tumor_fq1=tumor_fq1,
            tumor_fq2=tumor_fq2,

            normal_fq1=None if normal_fq1.lower() == 'none' else normal_fq1,
            normal_fq2=None if normal_fq2.lower() == 'none' else normal_fq2,

            clip_r1_5_prime=clip_r1_5_prime,
            clip_r2_5_prime=clip_r2_5_prime,
            read_aligner=read_aligner,
            skip_mark_duplicates=skip_mark_duplicates,
            bqsr_known_variant_vcf=None if bqsr_known_variant_vcf.lower() == 'none' else bqsr_known_variant_vcf,
            discard_bam=discard_bam,

            variant_callers=variant_callers.split(','),
            skip_variant_calling=skip_variant_calling,
            call_region_bed=None if call_region_bed.lower() == 'none' else call_region_bed,
            panel_of_normal_vcf=None if panel_of_normal_vcf.lower() == 'none' else panel_of_normal_vcf,
            germline_resource_vcf=None if germline_resource_vcf.lower() == 'none' else germline_resource_vcf,

            variant_flagging_criteria=None if variant_flagging_criteria.lower() == 'none' else variant_flagging_criteria,
            variant_removal_flags=[] if variant_removal_flags.lower() == 'none' else variant_removal_flags.split(','),

            min_snv_callers=min_snv_callers,
            min_indel_callers=min_indel_callers,

            annotator=annotator,
            skip_variant_annotation=skip_variant_annotation,
            vep_db_tar_gz=None if vep_db_tar_gz.lower() == 'none' else vep_db_tar_gz,
            vep_db_type=vep_db_type,
            vep_buffer_size=vep_buffer_size,
            dbnsfp_resource=None if dbnsfp_resource.lower() == 'none' else dbnsfp_resource,
            cadd_resource=None if cadd_resource.lower() == 'none' else cadd_resource,
            clinvar_vcf_gz=None if clinvar_vcf_gz.lower() == 'none' else clinvar_vcf_gz,
            dbsnp_vcf_gz=None if dbsnp_vcf_gz.lower() == 'none' else dbsnp_vcf_gz,
            snpsift_dbnsfp_txt_gz=None if snpsift_dbnsfp_txt_gz.lower() == 'none' else snpsift_dbnsfp_txt_gz,

            skip_cnv=skip_cnv,
            exome_target_bed=None if exome_target_bed.lower() == 'none' else exome_target_bed,
            cnvkit_annotate_txt=None if cnvkit_annotate_txt.lower() == 'none' else cnvkit_annotate_txt
        )

    def annotate(
            self,
            vcf: str,
            ref_fa: str,

            outdir: str,
            threads: str,
            debug: bool,

            annotator: str,
            vep_db_tar_gz: str,
            vep_db_type: str,
            vep_buffer_size: int,
            dbnsfp_resource: str,
            cadd_resource: str,
            clinvar_vcf_gz: str,
            dbsnp_vcf_gz: str,
            snpsift_dbnsfp_txt_gz: str):

        self.config_settings(outdir=outdir, threads=int(threads), debug=debug)

        Annotation(self.settings).main(
            vcf=vcf,
            ref_fa=ref_fa,

            annotator=annotator,
            vep_db_tar_gz=None if vep_db_tar_gz.lower() == 'none' else vep_db_tar_gz,
            vep_db_type=vep_db_type,
            vep_buffer_size=vep_buffer_size,
            dbnsfp_resource=None if dbnsfp_resource.lower() == 'none' else dbnsfp_resource,
            cadd_resource=None if cadd_resource.lower() == 'none' else cadd_resource,
            clinvar_vcf_gz=None if clinvar_vcf_gz.lower() == 'none' else clinvar_vcf_gz,
            dbsnp_vcf_gz=None if dbsnp_vcf_gz.lower() == 'none' else dbsnp_vcf_gz,
            snpsift_dbnsfp_txt_gz=None if snpsift_dbnsfp_txt_gz.lower() == 'none' else snpsift_dbnsfp_txt_gz,
        )

    def config_settings(self, outdir: str, threads: int, debug: bool):
        self.settings = Settings(
            workdir=get_temp_path(prefix='./somatic_pipeline_workdir_'),
            outdir=outdir,
            threads=threads,
            debug=debug,
            mock=False)

        for d in [self.settings.workdir, self.settings.outdir]:
            os.makedirs(d, exist_ok=True)
