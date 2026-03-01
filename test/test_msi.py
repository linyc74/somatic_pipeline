import shutil
from os.path import exists
from somatic_pipeline.msi import MSIsensor, MANTIS
from .setup import TestCase


class TestMSIsensor(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        for file in ['chr9.fa', 'tumor-sorted.bam', 'normal-sorted.bam']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')
        
        MSIsensor(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.workdir}/tumor-sorted.bam',
            normal_bam=f'{self.workdir}/normal-sorted.bam',
            bed_file=f'{self.indir}/chr9-exome-probes.bed'
        )

        for suffix in ['', '_dis', '_germline', '_somatic']:
            self.assertTrue(exists(f'{self.outdir}/msi/msisensor{suffix}'))

    def test_without_bed_file(self):
        for file in ['chr9.fa', 'tumor-sorted.bam', 'normal-sorted.bam']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')
        
        MSIsensor(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.workdir}/tumor-sorted.bam',
            normal_bam=f'{self.workdir}/normal-sorted.bam',
            bed_file=None
        )

        for suffix in ['', '_dis', '_germline', '_somatic']:
            self.assertTrue(exists(f'{self.outdir}/msi/msisensor{suffix}'))


class TestMANTIS(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):   
        for file in ['chr9.fa', 'tumor-sorted.bam', 'normal-sorted.bam']:
            shutil.copy(f'{self.indir}/{file}', f'{self.workdir}/{file}')

        MANTIS(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.workdir}/tumor-sorted.bam',
            normal_bam=f'{self.workdir}/normal-sorted.bam',
        )
