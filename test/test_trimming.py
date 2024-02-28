from somatic_pipeline.trimming import Trimming, RemoveUmiAndAdapter, strip_mate_3prime_umi
from .setup import TestCase


class TestTrimming(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_no_umi(self):
        tumor_fq1, tumor_fq2, normal_fq1, normal_fq2 = Trimming(self.settings).main(
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz',
            umi_length=0,
            clip_r1_5_prime=7,
            clip_r2_5_prime=7,
        )
        self.assertFileExists(f'{self.workdir}/tumor.1_val_1.fq.gz', tumor_fq1)
        self.assertFileExists(f'{self.workdir}/tumor.2_val_2.fq.gz', tumor_fq2)
        self.assertFileExists(f'{self.workdir}/normal.1_val_1.fq.gz', normal_fq1)
        self.assertFileExists(f'{self.workdir}/normal.2_val_2.fq.gz', normal_fq2)

    def test_trim_umi(self):
        tumor_fq1, tumor_fq2, normal_fq1, normal_fq2 = Trimming(self.settings).main(
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz',
            umi_length=7,
            clip_r1_5_prime=0,
            clip_r2_5_prime=0,
        )
        self.assertFileExists(f'{self.workdir}/tumor.1_umi_adapter_removed_val_1.fq.gz', tumor_fq1)
        self.assertFileExists(f'{self.workdir}/tumor.2_umi_adapter_removed_val_2.fq.gz', tumor_fq2)
        self.assertFileExists(f'{self.workdir}/normal.1_umi_adapter_removed_val_1.fq.gz', normal_fq1)
        self.assertFileExists(f'{self.workdir}/normal.2_umi_adapter_removed_val_2.fq.gz', normal_fq2)


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
            gz=True
        )
        self.assertFileExists(f'{self.workdir}/tumor.1_umi_adapter_removed.fastq.gz', fq1)
        self.assertFileExists(f'{self.workdir}/tumor.2_umi_adapter_removed.fastq.gz', fq2)


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

