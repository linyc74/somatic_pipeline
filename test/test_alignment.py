from somatic_pipeline.alignment import BwaIndexer, BwaAligner, Bowtie2Indexer, Bowtie2Aligner
from .setup import TestCase


class TestBWA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        index = BwaIndexer(self.settings).main(
            fna=f'{self.indir}/chr9.fa'
        )
        actual = BwaAligner(self.settings).main(
            index=index,
            fq1=f'{self.indir}/tumor.1.fq.gz',
            fq2=f'{self.indir}/tumor.2.fq.gz',
            sample_name='tumor',
            discard_bam=False
        )
        expected = f'{self.outdir}/tumor-sorted.bam'
        self.assertFileExists(expected, actual)


class TestBowtie2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        index = Bowtie2Indexer(self.settings).main(
            fna=f'{self.indir}/chr9.fa'
        )
        actual = Bowtie2Aligner(self.settings).main(
            index=index,
            fq1=f'{self.indir}/tumor.1.fq.gz',
            fq2=f'{self.indir}/tumor.2.fq.gz',
            sample_name='tumor',
            discard_bam=True
        )
        expected = f'{self.workdir}/tumor-sorted.bam'
        self.assertFileExists(expected, actual)
