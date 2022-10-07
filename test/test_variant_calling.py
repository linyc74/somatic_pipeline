import shutil
from somatic_pipeline.variant_calling import VariantCalling
from .setup import TestCase


class TestVariantCalling(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        for file in ['chr9.fa', 'tumor-sorted.bam', 'normal-sorted.bam']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')
        self.ref_fa = f'{self.workdir}/chr9.fa'
        self.tumor_bam = f'{self.workdir}/tumor-sorted.bam'
        self.normal_bam = f'{self.workdir}/normal-sorted.bam'

    def tearDown(self):
        self.tear_down()

    def test_mutect2_tn_paired(self):
        actual = VariantCalling(self.settings).main(
            variant_caller='mutect2',
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            filter_variants=True,
            variant_removal_flags=['orientation'],
        )
        expected = f'{self.workdir}/raw-filter-mutect-calls-variant-removal.vcf'
        self.assertFileExists(expected, actual)

    def test_mutect2_tumor_only(self):
        actual = VariantCalling(self.settings).main(
            variant_caller='mutect2',
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=None,
            panel_of_normal_vcf=f'{self.indir}/22_0830_combine_pon_chr9.vcf.gz',
            germline_resource_vcf=f'{self.indir}/af-only-gnomad.hg38.chr9.vcf.gz',
            filter_variants=True,
            variant_removal_flags=['orientation', 'panel_of_normals'],
        )
        expected = f'{self.workdir}/raw-filter-mutect-calls-variant-removal.vcf'
        self.assertFileExists(expected, actual)

    def test_haplotype_caller(self):
        actual = VariantCalling(self.settings).main(
            variant_caller='haplotype-caller',
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=None,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            filter_variants=True,
            variant_removal_flags=[]
        )
        expected = f'{self.workdir}/raw-snp-indel-flagged.vcf'
        self.assertFileExists(expected, actual)

    def test_muse(self):
        actual = VariantCalling(self.settings).main(
            variant_caller='muse',
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            filter_variants=False,
            variant_removal_flags=[],
        )
        expected = f'{self.workdir}/raw.vcf'
        self.assertFileExists(expected, actual)

    def test_varscan(self):
        actual = VariantCalling(self.settings).main(
            variant_caller='varscan',
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            filter_variants=False,
            variant_removal_flags=[],
        )
        expected = f'{self.workdir}/raw.vcf'
        self.assertFileExists(expected, actual)

    def test_raise_assertion_error(self):
        with self.assertRaises(AssertionError):
            VariantCalling(self.settings).main(
                variant_caller='invalid-caller',
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                panel_of_normal_vcf=None,
                germline_resource_vcf=None,
                filter_variants=False,
                variant_removal_flags=[]
            )

        for caller in ['muse', 'varscan']:
            with self.assertRaises(AssertionError):
                VariantCalling(self.settings).main(
                    variant_caller=caller,
                    ref_fa=self.ref_fa,
                    tumor_bam=self.tumor_bam,
                    normal_bam=None,
                    panel_of_normal_vcf=None,
                    germline_resource_vcf=None,
                    filter_variants=False,
                    variant_removal_flags=[]
                )
