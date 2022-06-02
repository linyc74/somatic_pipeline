from somatic_pipeline.annotation import Annotation, SnpEff, SnpSiftDbNSFP, SnpSiftAnnotate
from .setup import TestCase


class TestAnnotation(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = Annotation(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,
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


class TestSnpSiftAnnotate(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SnpSiftAnnotate(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            db_vcf_gz=f'{self.indir}/clinvar.vcf.gz'
        )
        expected = f'{self.workdir}/raw-clinvar.vcf'
        self.assertFileExists(expected, actual)


class TestSnpSiftDbNSFP(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SnpSiftDbNSFP(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            dbnsfp_txt_gz=f'{self.indir}/22_0414_dbNSFP_chr9.txt.gz'
        )
        expected = f'{self.workdir}/raw-dbnsfp.vcf'
        self.assertFileExists(expected, actual)
