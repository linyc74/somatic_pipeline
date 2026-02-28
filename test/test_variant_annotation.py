import shutil
from somatic_pipeline.variant_annotation import VariantAnnotation, VEP, SnpSiftAnnotate, AssertDbnsfpResourceFilenameForVEP
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
        actual = VariantAnnotation(self.settings).main(
            vcf=f'{self.indir}/picked-variants.vcf',
            ref_fa=self.ref_fa,
            clinvar_vcf_gz=f'{self.indir}/clinvar.vcf.gz',
            dbsnp_vcf_gz=None,
            vep_db_tar_gz=f'{self.indir}/homo_sapiens_merged_vep_113_GRCh38_chr9.tar.gz',
            vep_db_type='merged',
            vep_buffer_size=100,
            cadd_resource=None,
            dbnsfp_resource=f'{self.indir}/22_0414_dbNSFP_chr9_4.1a.txt.gz',
        )
        expected = f'{self.workdir}/picked-variants-vep-clinvar.vcf'
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
            vcf=f'{self.indir}/picked-variants.vcf',
            ref_fa=self.ref_fa,
            vep_db_tar_gz=f'{self.indir}/homo_sapiens_merged_vep_113_GRCh38_chr9.tar.gz',
            vep_db_type='merged',
            vep_buffer_size=100,
            cadd_resource=None,
            dbnsfp_resource=f'{self.indir}/22_0414_dbNSFP_chr9_4.1a.txt.gz',
        )
        expected = f'{self.workdir}/picked-variants-vep.vcf'
        self.assertEqual(expected, actual)
        self.assertFileExists(expected, actual)


class TestSnpSiftAnnotate(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SnpSiftAnnotate(self.settings).main(
            vcf=f'{self.indir}/picked-variants.vcf',
            resource_vcf_gz=f'{self.indir}/clinvar.vcf.gz'
        )
        expected = f'{self.workdir}/picked-variants-clinvar.vcf'
        self.assertFileExists(expected, actual)


class TestAssertDbnsfpResourceFilenameForVEP(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_right_filename(self):
        # contain '3.' or '4.'
        AssertDbnsfpResourceFilenameForVEP(self.settings).main(fpath='dbNSFP3._a.gz')
        AssertDbnsfpResourceFilenameForVEP(self.settings).main(fpath='dbNSFP4._a.gz')

    def test_wrong_filename(self):
        with self.assertRaises(AssertionError):
            # not contain '3.' or '4.'
            AssertDbnsfpResourceFilenameForVEP(self.settings).main(fpath='dbNSFP_.1a.gz')
