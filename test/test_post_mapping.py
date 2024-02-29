import shutil
from os.path import exists
from somatic_pipeline.post_mapping import BQSR
from somatic_pipeline.post_mapping import SamtoolsStats
from somatic_pipeline.post_mapping import GATKMarkDuplicates
from .setup import TestCase


class TestGATKMarkDuplicates(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = GATKMarkDuplicates(self.settings).main(
            bam=f'{self.indir}/tumor-sorted.bam'
        )
        expected = f'{self.workdir}/tumor-sorted-mark-duplicates.bam'
        self.assertFileExists(expected, actual)


class TestBQSR(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        shutil.copy(f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa')
        self.ref_fa = f'{self.workdir}/chr9.fa'
        self.tumor_bam = f'{self.indir}/tumor-sorted.bam'
        self.normal_bam = f'{self.indir}/normal-sorted.bam'

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        BQSR(self.settings).main(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            ref_fa=self.ref_fa,
            known_variant_vcf=f'{self.indir}/known-variants.vcf'
        )


class TestSamtoolsStats(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        SamtoolsStats(self.settings).main(
            bam=f'{self.indir}/tumor-sorted-mark-duplicates.bam')
        file = f'{self.outdir}/mapping-stats/tumor-sorted-mark-duplicates.txt'
        self.assertTrue(exists(file))
