import shutil
from somatic_pipeline.cnv import CNV
from .setup import TestCase


class TestCNV(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        for f in ['chr9.fa', 'tumor-sorted.bam', 'normal-sorted.bam']:
            shutil.copy(f'{self.indir}/{f}', f'{self.workdir}/{f}')

        CNV(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.workdir}/tumor-sorted.bam',
            normal_bam=f'{self.workdir}/normal-sorted.bam',
            exome_target_bed=f'{self.indir}/chr9-exome-probes.bed',
            annotate_txt=f'{self.indir}/refFlat.txt',
            segmentation_threshold=1e-6,
        )
