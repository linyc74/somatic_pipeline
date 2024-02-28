from somatic_pipeline.trimming import Trimming
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
