from somatic_pipeline.trimming import Trimming, RemoveUmiAndAdapter, strip_mate_3prime_umi
from .setup import TestCase


class TestTrimming(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        tumor_fq1, tumor_fq2, normal_fq1, normal_fq2 = Trimming(self.settings).main(
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=None,
            normal_fq2=None,
            clip_r1_5_prime=1,
            clip_r2_5_prime=1,
        )
        self.assertFileExists(f'{self.workdir}/tumor.1_val_1.fq.gz', tumor_fq1)
        self.assertFileExists(f'{self.workdir}/tumor.2_val_2.fq.gz', tumor_fq2)
        self.assertIsNone(normal_fq1)
        self.assertIsNone(normal_fq2)


class TestRemoveUmiAndAdapter(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        fq1, fq2 = RemoveUmiAndAdapter(self.settings).main(
            fq1=f'{self.indir}/tumor.1.fq.gz',
            fq2=f'{self.indir}/tumor.2.fq.gz',
            umi_length=7,
            gz=False
        )
        self.assertFileExists(f'{self.workdir}/umi_adapter_removed_R1.fastq', fq1)
        self.assertFileExists(f'{self.workdir}/umi_adapter_removed_R2.fastq', fq2)


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

