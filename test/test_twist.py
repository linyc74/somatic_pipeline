from somatic_pipeline.twist import TwistUMIConsensusMapping, RemoveTwistUMI, strip_mate_3prime_umi
from .setup import TestCase


class TestTwistUMIConsensusMapping(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        TwistUMIConsensusMapping(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa',
            fq1=f'{self.indir}/twist_tiny_tumor_R1.fastq.gz',
            fq2=f'{self.indir}/twist_tiny_tumor_R2.fastq.gz',
            sample_name='twist_tiny_tumor'
        )


class TestRemoveTwistUMI(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        fq1, fq2 = RemoveTwistUMI(self.settings).main(
            fq1=f'{self.indir}/twist_tiny_tumor_R1.fastq.gz',
            fq2=f'{self.indir}/twist_tiny_tumor_R2.fastq.gz',
            gz=False
        )
        self.assertFileExists(f'{self.workdir}/twist_umi_removed_R1.fastq', fq1)
        self.assertFileExists(f'{self.workdir}/twist_umi_removed_R2.fastq', fq2)


class TestFunctions(TestCase):

    def test_has_umi(self):
        actual = strip_mate_3prime_umi(
            read='A' * 20,
            mate='T' * 20 + 'ACGTACGT'
        )
        self.assertEqual('T' * 20, actual)

    def test_no_umi(self):
        actual = strip_mate_3prime_umi(
            read='A' * 20,
            mate='T' * 20
        )
        self.assertEqual('T' * 20, actual)

    def test_no_overlap(self):
        actual = strip_mate_3prime_umi(
            read='A' * 20,
            mate='G' * 20
        )
        self.assertEqual('G' * 20, actual)
