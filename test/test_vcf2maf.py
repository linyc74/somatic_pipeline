import shutil
from somatic_pipeline.vcf2maf import Vcf2Maf
from .setup import TestCase


class TestVcf2Maf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        shutil.copy(f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa')
        self.ref_fa = f'{self.workdir}/chr9.fa'
        self.annotated_vcf = f'{self.indir}/annotated.vcf'

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = Vcf2Maf(self.settings).main(
            annotated_vcf=self.annotated_vcf,
            ref_fa=self.ref_fa,
            variant_caller='mutect2',
        )
        expected = f'{self.indir}/variants.maf'
        self.assertFileEqual(expected, actual)
