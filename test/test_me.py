from .setup import TestCase
from gatk_pipeline.trimming import TrimGalore


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_me(self):
        fq1, fq2 = TrimGalore(self.settings).main(
            fq1=f'{self.indir}/normal.1.fq.gz',
            fq2=f'{self.indir}/normal.2.fq.gz'
        )
        self.assertFileExists(
            expected=f'{self.workdir}/normal.1_val_1.fq.gz',
            actual=fq1
        )
        self.assertFileExists(
            expected=f'{self.workdir}/normal.2_val_2.fq.gz',
            actual=fq2
        )
