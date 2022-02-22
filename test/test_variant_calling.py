import shutil
from somatic_pipeline.variant_calling import Mutect2TumorNormalPaired, Mutect2TumorOnly
from .setup import TestCase


class TestMutect2TumorNormalPaired(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        shutil.copy(
            f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa'  # ref_fa should be in workdir during runtime
        )
        actual = Mutect2TumorNormalPaired(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.indir}/tumor-sorted.bam',
            normal_bam=f'{self.indir}/normal-sorted.bam'
        )
        expected = f'{self.outdir}/raw.vcf'
        self.assertFileExists(expected, actual)


class TestMutect2TumorOnly(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        shutil.copy(
            f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa'  # ref_fa should be in workdir during runtime
        )
        actual = Mutect2TumorOnly(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.indir}/tumor-sorted.bam',
        )
        expected = f'{self.outdir}/raw.vcf'
        self.assertFileExists(expected, actual)
