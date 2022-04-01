import shutil
from somatic_pipeline.bqsr import RunBQSR
from .setup import TestCase


class TestRunBQSR(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        shutil.copy(f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa')
        self.ref_fa = f'{self.workdir}/chr9.fa'
        self.tumor_bam = f'{self.indir}/tumor-sorted.bam'

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        RunBQSR(self.settings).main(
            bam=self.tumor_bam,
            ref_fa=self.ref_fa,
            known_variant_vcf=f'{self.indir}/known-variants.vcf'
        )
