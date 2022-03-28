from somatic_pipeline.mark_duplicates import GATKMarkDuplicates
from .setup import TestCase


class TestGATKMarkDuplicates(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = GATKMarkDuplicates(self.settings).main(
            bam=f'{self.indir}/tumor-sorted.bam'
        )
        expected = f'{self.workdir}/tumor-sorted-mark-duplicates.bam'
        self.assertFileExists(expected, actual)
