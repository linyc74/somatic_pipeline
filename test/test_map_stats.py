from os.path import exists
from somatic_pipeline.map_stats import SamtoolsStats
from .setup import TestCase


class TestSamtoolsStats(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        SamtoolsStats(self.settings).main(
            bam=f'{self.indir}/tumor-sorted-mark-duplicates.bam')
        file = f'{self.outdir}/mapping-stats/tumor-sorted-mark-duplicates.txt'
        self.assertTrue(exists(file))
