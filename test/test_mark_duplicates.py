from somatic_pipeline.mark_duplicates import MarkDuplicates
from .setup import TestCase


class TestMarkDuplicates(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = MarkDuplicates(self.settings).main(
            bam_in=f'{self.indir}/tumor-sorted.bam'
        )
        expected = f'{self.outdir}/tumor-sorted-mark-duplicates.bam'
        self.assertFileExists(expected, actual)
