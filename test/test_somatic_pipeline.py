from .setup import TestCase
from somatic_pipeline.somatic_pipeline import SomaticPipeline


class TestSomaticPipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_tumor_normal_paired(self):
        SomaticPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa.gz',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',

            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz',

            clip_r1_5_prime=1,
            clip_r2_5_prime=1,
            read_aligner='bwa',
            skip_mark_duplicates=False,
            bqsr_known_variant_vcf=f'{self.indir}/known-variants.vcf',
            discard_bam=True,

            variant_callers=['mutect2', 'lofreq'],
            skip_variant_calling=False,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            vardict_call_region_bed=None,

            variant_removal_flags=[],

            min_snv_callers=1,
            min_indel_callers=1,

            annotator='vep',
            skip_variant_annotation=False,
            vep_db_tar_gz=f'{self.indir}/homo_sapiens_merged_vep_106_GRCh38_chr9.tar.gz',
            vep_db_type='merged',
            vep_buffer_size=100,
            dbnsfp_resource=f'{self.indir}/22_0414_dbNSFP_chr9_4.1a.txt.gz',
            cadd_resource=None,
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,
            snpsift_dbnsfp_txt_gz=f'{self.indir}/22_0414_dbNSFP_chr9_4.1a.txt.gz',

            skip_cnv=False,
            exome_target_bed=f'{self.indir}/chr9-exome-probes.bed',
            cnvkit_annotate_txt=f'{self.indir}/chr9-refFlat.txt'
        )

    def test_tumor_only(self):
        SomaticPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa.gz',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',

            normal_fq1=None,
            normal_fq2=None,

            clip_r1_5_prime=1,
            clip_r2_5_prime=1,
            read_aligner='bwa',
            skip_mark_duplicates=False,
            bqsr_known_variant_vcf=f'{self.indir}/known-variants.vcf',
            discard_bam=True,

            variant_callers=['mutect2', 'lofreq'],
            skip_variant_calling=False,
            panel_of_normal_vcf=f'{self.indir}/22_0830_combine_pon_chr9.vcf.gz',
            germline_resource_vcf=f'{self.indir}/af-only-gnomad.hg38.chr9.vcf.gz',
            vardict_call_region_bed=f'{self.indir}/chr9-exome-probes.bed',

            variant_removal_flags=['panel_of_normals', 'map_qual'],

            min_snv_callers=1,
            min_indel_callers=1,

            annotator='snpeff',
            skip_variant_annotation=False,
            vep_db_tar_gz=None,
            vep_db_type='merged',
            vep_buffer_size=100,
            dbnsfp_resource=None,
            cadd_resource=None,
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,
            snpsift_dbnsfp_txt_gz=None,

            skip_cnv=False,
            exome_target_bed=None,
            cnvkit_annotate_txt=None
        )
