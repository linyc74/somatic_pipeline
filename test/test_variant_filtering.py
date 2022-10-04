import shutil
from somatic_pipeline.variant_filtering import VariantFiltering
from .setup import TestCase


class TestVariantFiltering(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        for file in ['chr9.fa']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')
        self.ref_fa = f'{self.workdir}/chr9.fa'

    def tearDown(self):
        self.tear_down()

    def test_mutect2_filtering(self):
        actual = VariantFiltering(self.settings).main(
            vcf=f'{self.indir}/raw-mutect2.vcf',
            ref_fa=self.ref_fa,
            variant_caller='mutect2',
            variant_removal_flags=['panel_of_normals', 'map_qual']
        )
        expected = f'{self.workdir}/raw-mutect2-filter-mutect-calls-variant-removal.vcf'
        self.assertFileExists(expected, actual)

    def test_haplotype_filtering(self):
        actual = VariantFiltering(self.settings).main(
            vcf=f'{self.indir}/raw-haplotype.vcf',
            ref_fa=self.ref_fa,
            variant_caller='haplotype-caller',
            variant_removal_flags=['QD2']
        )
        expected = f'{self.workdir}/raw-haplotype-snp-indel-flagged-variant-removal.vcf'
        self.assertFileExists(expected, actual)

    def test_skip_for_invalid_caller(self):
        actual = VariantFiltering(self.settings).main(
            vcf=f'{self.indir}/raw-mutect2.vcf',
            ref_fa=self.ref_fa,
            variant_caller='muse',
            variant_removal_flags=[]
        )
        expected = f'{self.indir}/raw-mutect2.vcf'
        self.assertFileExists(expected, actual)
