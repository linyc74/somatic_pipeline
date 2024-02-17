from somatic_pipeline.mapping import Mapping
from .setup import TestCase


class TestMapping(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_bwa(self):
        tumor_bam, normal_bam = Mapping(self.settings).main(
            read_aligner='bwa',
            ref_fa=f'{self.indir}/chr9.fa',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz',
            discard_bam=False)
        self.assertFileExists(f'{self.outdir}/tumor-sorted.bam', tumor_bam)
        self.assertFileExists(f'{self.outdir}/normal-sorted.bam', normal_bam)

    def test_bowtie2(self):
        tumor_bam, normal_bam = Mapping(self.settings).main(
            read_aligner='bowtie2',
            ref_fa=f'{self.indir}/chr9.fa',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=None,
            normal_fq2=None,
            discard_bam=False)
        self.assertFileExists(f'{self.outdir}/tumor-sorted.bam', tumor_bam)
        self.assertIsNone(normal_bam)
