from somatic_pipeline.annotation import Annotation, SnpEff, SnpSift
from .setup import TestCase


class TestAnnotation(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = Annotation(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            snpsift_dbnsfp_txt_gz=f'{self.indir}/22_0414_dbNSFP_chr9.txt.gz',
        )
        expected = f'{self.outdir}/annotated.vcf'
        self.assertFileExists(expected, actual)


class TestSnpEff(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SnpEff(self.settings).main(
            vcf=f'{self.indir}/raw.vcf'
        )
        expected = f'{self.workdir}/raw-snpeff.vcf'
        self.assertFileExists(expected, actual)


class TestSnpSift(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SnpSift(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            dbnsfp_txt_gz=f'{self.indir}/22_0414_dbNSFP_chr9.txt.gz'
        )
        expected = f'{self.workdir}/raw-snpsift.vcf'
        self.assertFileExists(expected, actual)
