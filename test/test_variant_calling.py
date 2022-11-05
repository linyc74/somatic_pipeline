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

    # def tearDown(self):
    #     self.tear_down()

    def test_tn_paired(self):
        variant_callers = [
            'lofreq',
            'somatic-sniper',
            'vardict',
            'varscan',
            'muse',
            'mutect2',
        ]
        actual = VariantCalling(self.settings).main(
            variant_callers=variant_callers,
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            vardict_call_region_bed=None,
            variant_removal_flags=[],
        )
        expected = [
            f'{self.outdir}/callers/lofreq.vcf',
            f'{self.outdir}/callers/somatic-sniper.vcf',
            f'{self.outdir}/callers/vardict.vcf',
            f'{self.outdir}/callers/varscan.vcf',
            f'{self.outdir}/callers/muse.vcf',
            f'{self.outdir}/callers/mutect2.vcf',
        ]
        for a, e in zip(actual, expected):
            self.assertFileExists(e, a)

    def test_tumor_only(self):
        variant_callers = [
            'lofreq',
            'vardict',
            'haplotype-caller',
            'mutect2',
        ]
        actual = VariantCalling(self.settings).main(
            variant_callers=variant_callers,
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=None,
            germline_resource_vcf=None,
            vardict_call_region_bed=None,
            variant_removal_flags=[],
        )
        expected = [
            f'{self.outdir}/callers/lofreq.vcf',
            f'{self.outdir}/callers/vardict.vcf',
            f'{self.outdir}/callers/haplotype-caller.vcf',
            f'{self.outdir}/callers/mutect2.vcf',
        ]
        for a, e in zip(actual, expected):
            self.assertFileExists(e, a)

    def test_raise_assertion_error(self):
        with self.assertRaises(AssertionError):
            VariantCalling(self.settings).main(
                variant_callers=['invalid-caller'],
                ref_fa=self.ref_fa,
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                panel_of_normal_vcf=None,
                germline_resource_vcf=None,
                vardict_call_region_bed=None,
                variant_removal_flags=[]
            )
