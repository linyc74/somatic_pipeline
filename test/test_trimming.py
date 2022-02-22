from somatic_pipeline.trimming import TrimGalore
from .setup import TestCase


class TestTrimGalore(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        trimmed_fq1, trimmed_fq2 = TrimGalore(self.settings).main(
            fq1=f'{self.indir}/tumor.1.fq.gz',
            fq2=f'{self.indir}/tumor.2.fq.gz',
        )
        for expected, actual in [
            (f'{self.workdir}/tumor.1_val_1.fq.gz', trimmed_fq1),
            (f'{self.workdir}/tumor.2_val_2.fq.gz', trimmed_fq2),
        ]:
            self.assertFileExists(expected, actual)

