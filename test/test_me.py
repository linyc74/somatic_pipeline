from .setup import TestCase
from gatk_pipeline.mapping import BwaIndex, BwaMem


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_me(self):
        index = BwaIndex(self.settings).main(
            fna=f'{self.indir}/chr9.fa')

        sorted_bam = BwaMem(self.settings).main(
            index=index,
            fq1=f'{self.indir}/normal.1.fq.gz',
            fq2=f'{self.indir}/normal.2.fq.gz',
            sample_name='normal')

        self.assertFileExists(
            expected=f'{self.outdir}/normal_sorted.bam',
            actual=sorted_bam
        )
