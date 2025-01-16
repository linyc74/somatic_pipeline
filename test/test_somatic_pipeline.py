from .setup import TestCase
from somatic_pipeline.somatic_pipeline import SomaticPipeline


class TestSomaticPipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_tumor_normal_paired(self):
        SomaticPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa.gz',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',

            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz',

            umi_length=3,
            clip_r1_5_prime=1,
            clip_r2_5_prime=1,
            read_aligner='bwa',
            skip_mark_duplicates=False,
            bqsr_known_variant_vcf=f'{self.indir}/known-variants.vcf',
            discard_bam=True,

            variant_callers=['mutect2', 'vardict', 'lofreq', 'muse', 'varscan', 'somatic-sniper'],
            skip_variant_calling=False,
            panel_of_normal_vcf=f'{self.indir}/22_0830_combine_pon_chr9.vcf.gz',
            germline_resource_vcf=f'{self.indir}/af-only-gnomad.hg38.chr9.vcf.gz',
            call_region_bed=f'{self.indir}/chr9-exome-probes.bed',

            variant_flagging_criteria='low_depth:DP<10',
            variant_removal_flags=[],
            only_pass=True,

            min_snv_callers=1,
            min_indel_callers=1,

            skip_variant_annotation=False,
            vep_db_tar_gz=f'{self.indir}/homo_sapiens_merged_vep_106_GRCh38_chr9.tar.gz',
            vep_db_type='merged',
            vep_buffer_size=100,
            dbnsfp_resource=f'{self.indir}/22_0414_dbNSFP_chr9_4.1a.txt.gz',
            cadd_resource=None,
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,

            pcgr_ref_data_tgz=f'{self.indir}/homo_sapiens_vep_112_GRCh38_chr9_chr22.tar.gz',
            pcgr_vep_tar_gz=f'{self.indir}/pcgr_ref_data.20240927.grch38.tgz',
            pcgr_tumor_site=12,
            pcgr_tmb_target_size_mb=34,
            pcgr_tmb_display='coding_and_silent',
        )

    def test_tumor_only(self):
        SomaticPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa.gz',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',

            normal_fq1=None,
            normal_fq2=None,

            umi_length=0,
            clip_r1_5_prime=1,
            clip_r2_5_prime=1,
            read_aligner='bwa',
            skip_mark_duplicates=False,
            bqsr_known_variant_vcf=f'{self.indir}/known-variants.vcf',
            discard_bam=True,

            variant_callers=['mutect2', 'haplotype-caller', 'vardict', 'lofreq'],
            skip_variant_calling=False,
            panel_of_normal_vcf=f'{self.indir}/22_0830_combine_pon_chr9.vcf.gz',
            germline_resource_vcf=f'{self.indir}/af-only-gnomad.hg38.chr9.vcf.gz',
            call_region_bed=f'{self.indir}/chr9-exome-probes.bed',

            variant_flagging_criteria='low_depth:DP<10',
            variant_removal_flags=['panel_of_normals', 'map_qual'],
            only_pass=True,

            min_snv_callers=1,
            min_indel_callers=1,

            skip_variant_annotation=False,
            vep_db_tar_gz=None,
            vep_db_type='merged',
            vep_buffer_size=100,
            dbnsfp_resource=None,
            cadd_resource=None,
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,

            pcgr_ref_data_tgz=None,
            pcgr_vep_tar_gz=None,
            pcgr_tumor_site=12,
            pcgr_tmb_target_size_mb=34,
            pcgr_tmb_display='coding_and_silent',
        )
