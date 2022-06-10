import shutil
from somatic_pipeline.annotation import Annotation, SnpEff, SnpSiftDbNSFP, SnpSiftAnnotate, VEP
from .setup import TestCase


class TestAnnotation(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def tearDown(self):
        self.tear_down()

    def copy_and_set_files(self):
        for file in ['chr9.fa']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')
        self.ref_fa = f'{self.workdir}/chr9.fa'

    def test_main(self):
        actual = Annotation(self.settings).main(
            annotator='vep',
            vcf=f'{self.indir}/raw.vcf',
            ref_fa=self.ref_fa,
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,
            snpsift_dbnsfp_txt_gz=f'{self.indir}/22_0414_dbNSFP_chr9.txt.gz',
            vep_db_tar_gz=f'{self.indir}/homo_sapiens_merged_vep_106_GRCh38_chr9.tar.gz',
            vep_db_type='merged'
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


class TestVEP(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def tearDown(self):
        self.tear_down()

    def copy_and_set_files(self):
        for file in ['chr9.fa']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')
        self.ref_fa = f'{self.workdir}/chr9.fa'

    def test_main(self):
        actual = VEP(self.settings).main(
            vcf=f'{self.indir}/raw.vcf',
            ref_fa=self.ref_fa,
            vep_db_tar_gz=f'{self.indir}/homo_sapiens_merged_vep_106_GRCh38_chr9.tar.gz',
            vep_db_type='merged'
        )
        expected = f'{self.workdir}/raw-vep.vcf'
        self.assertEqual(expected, actual)
        self.assertFileExists(expected, actual)
