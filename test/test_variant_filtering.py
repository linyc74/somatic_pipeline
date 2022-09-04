import shutil
from somatic_pipeline.variant_filtering import Mutect2VariantFiltering
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

    def test_main(self):
        actual = Mutect2VariantFiltering(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            ref_fa=self.ref_fa,
            variant_caller='mutect2',
            variant_removal_flags=['panel_of_normals', 'map_qual']
        )
        expected = f'{self.workdir}/raw-filter-mutect-calls-variant-removal.vcf'
        self.assertFileExists(expected, actual)

    def test_not_valid_caller(self):
        vcf = f'{self.indir}/raw.vcf'
        actual = Mutect2VariantFiltering(self.settings).main(
            vcf=vcf,
            ref_fa=self.ref_fa,
            variant_caller='muse',
            variant_removal_flags=[]
        )
        self.assertFileExists(vcf, actual)
