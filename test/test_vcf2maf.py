import shutil
from somatic_pipeline.vcf2maf import Vcf2Maf
from .setup import TestCase


class TestVcf2Maf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)
        self.copy_and_set_files()

    def copy_and_set_files(self):
        shutil.copy(f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa')
        self.ref_fa = f'{self.workdir}/chr9.fa'

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        for vcf in [
            'tn-paired-lofreq.vcf',
            'tn-paired-muse.vcf',
            'tn-paired-mutect2.vcf',
            'tn-paired-picked-variants.vcf',
            'tn-paired-somatic-sniper.vcf',
            'tn-paired-vardict.vcf',
            'tn-paired-varscan.vcf',
            'tumor-only-haplotype-caller.vcf',
            'tumor-only-lofreq.vcf',
            'tumor-only-mutect2.vcf',
            'tumor-only-picked-variants.vcf',
            'tumor-only-vardict.vcf',
        ]:
            with self.subTest(vcf=vcf):
                actual = Vcf2Maf(self.settings).main(
                    vcf=f'{self.indir}/{vcf}',
                    ref_fa=self.ref_fa,
                    dstdir=self.outdir
                )
                expected = f'{self.outdir}/{vcf[:-4]}.maf'
                self.assertFileExists(expected, actual)
