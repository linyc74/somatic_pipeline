import os
from typing import Optional
from .template import Settings
from .somatic_pipeline import SomaticPipeline


class Main:

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
    clinvar_vcf_gz: Optional[str]
    dbsnp_vcf_gz: Optional[str]
    snpsift_dbnsfp_txt_gz: Optional[str]
    discard_bam: bool
    skip_mark_duplicates: bool
    skip_variant_calling: bool
    skip_cnv: bool

    settings: Settings

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: str,
            normal_fq2: str,
            read_aligner: str,
            variant_caller: str,
            exome_target_bed: str,
            cnvkit_annotate_txt: str,
            panel_of_normal_vcf: str,
            bqsr_known_variant_vcf: str,
            clinvar_vcf_gz: str,
            dbsnp_vcf_gz: str,
            snpsift_dbnsfp_txt_gz: str,
            discard_bam: bool,
            skip_mark_duplicates: bool,
            skip_variant_calling: bool,
            skip_cnv: bool,
            outdir: str,
            threads: str,
            debug: bool):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = None if normal_fq1.lower() == 'none' else normal_fq1
        self.normal_fq2 = None if normal_fq2.lower() == 'none' else normal_fq2
        self.read_aligner = read_aligner
        self.variant_caller = variant_caller
        self.exome_target_bed = None if exome_target_bed.lower() == 'none' else exome_target_bed
        self.cnvkit_annotate_txt = None if cnvkit_annotate_txt.lower() == 'none' else cnvkit_annotate_txt
        self.panel_of_normal_vcf = None if panel_of_normal_vcf.lower() == 'none' else panel_of_normal_vcf
        self.bqsr_known_variant_vcf = None if bqsr_known_variant_vcf.lower() == 'none' else bqsr_known_variant_vcf
        self.clinvar_vcf_gz = None if clinvar_vcf_gz.lower() == 'none' else clinvar_vcf_gz
        self.dbsnp_vcf_gz = None if dbsnp_vcf_gz.lower() == 'none' else dbsnp_vcf_gz
        self.snpsift_dbnsfp_txt_gz = None if snpsift_dbnsfp_txt_gz.lower() == 'none' else snpsift_dbnsfp_txt_gz
        self.discard_bam = discard_bam
        self.skip_mark_duplicates = skip_mark_duplicates
        self.skip_variant_calling = skip_variant_calling
        self.skip_cnv = skip_cnv

        self.settings = Settings(
            workdir='./somatic_pipeline_workdir',
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
            variant_caller=self.variant_caller,
            exome_target_bed=self.exome_target_bed,
            cnvkit_annotate_txt=self.cnvkit_annotate_txt,
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            bqsr_known_variant_vcf=self.bqsr_known_variant_vcf,
            clinvar_vcf_gz=self.clinvar_vcf_gz,
            dbsnp_vcf_gz=self.dbsnp_vcf_gz,
            snpsift_dbnsfp_txt_gz=self.snpsift_dbnsfp_txt_gz,
            discard_bam=self.discard_bam,
            skip_mark_duplicates=self.skip_mark_duplicates,
            skip_variant_calling=self.skip_variant_calling,
            skip_cnv=self.skip_cnv)
