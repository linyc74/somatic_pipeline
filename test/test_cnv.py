import shutil
from somatic_pipeline.cnv import CNVkit
from .setup import TestCase


class TestCNVkit(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        for f in ['chr9.fa', 'tumor-sorted.bam', 'normal-sorted.bam']:
            shutil.copy(f'{self.indir}/{f}', f'{self.workdir}/{f}')
        self.ref_fa = f'{self.workdir}/chr9.fa'
        self.tumor_bam = f'{self.workdir}/tumor-sorted.bam'
        self.normal_bam = f'{self.workdir}/normal-sorted.bam'

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        CNVkit(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            gene_annotation_gff=f'{self.indir}/chr9.gff3',
            exome_target_bed=None
        )