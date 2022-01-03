from .setup import TestCase
from gatk_pipeline.gatk_pipeline import GATKPipeline


class TestGATKPipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        GATKPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz')
