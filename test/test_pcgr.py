from somatic_pipeline.pcgr import PCGR
from .setup import TestCase


class TestPCGR(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCGR(self.settings).main(
            vcf=f'{self.indir}/mutect2.vcf.gz',
            pcgr_ref_data_tgz=f'{self.indir}/pcgr_ref_data.20240927.grch38.tgz',
            pcgr_vep_tar_gz=f'{self.indir}/homo_sapiens_vep_112_GRCh38_chr9_chr22.tar.gz',
            vep_buffer_size=100,
            pcgr_tumor_site=12,
            pcgr_tmb_target_size_mb=34,
            pcgr_tmb_display='coding_and_silent'
        )
