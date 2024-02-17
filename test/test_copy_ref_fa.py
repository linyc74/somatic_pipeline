from os.path import exists
from somatic_pipeline.somatic_pipeline import CopyRefFa
from .setup import TestCase


class TestCopyRefFa(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = CopyRefFa(self.settings).main(ref_fa=f'{self.indir}/chr9.fa.gz')
        expected = f'{self.workdir}/chr9.fa'
        self.assertFileExists(expected, actual)
        for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
            with self.subTest(ext=ext):
                self.assertTrue(exists(f'{self.workdir}/bwa-index{ext}'))
