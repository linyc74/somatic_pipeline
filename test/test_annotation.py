from somatic_pipeline.annotation import SnpEff
from .setup import TestCase


class TestSnpEff(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SnpEff(self.settings).main(
            vcf=f'{self.indir}/raw.vcf'
        )
        expected = f'{self.outdir}/annotated.vcf'
        self.assertFileExists(expected, actual)
