import shutil
from somatic_pipeline.annotation import SnpEff
from somatic_pipeline.trimming import TrimGalore
from somatic_pipeline.somatic_pipeline import CopyRefFa
from somatic_pipeline.variant_calling import Mutect2TumorNormalPaired, Mutect2TumorOnly
from .setup import TestCase


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def __test_copy_ref_fa(self):
        actual = CopyRefFa(self.settings).main(ref_fa=f'{self.indir}/chr9.fa.gz')
        expected = f'{self.workdir}/chr9.fa'
        self.assertFileExists(expected, actual)

    def __test_trim_galore(self):
        trimmed_fq1, trimmed_fq2 = TrimGalore(self.settings).main(
            fq1=f'{self.indir}/tumor.1.fq.gz',
            fq2=f'{self.indir}/tumor.2.fq.gz',
        )
        for expected, actual in [
            (f'{self.workdir}/tumor.1_val_1.fq.gz', trimmed_fq1),
            (f'{self.workdir}/tumor.2_val_2.fq.gz', trimmed_fq2),
        ]:
            self.assertFileExists(expected, actual)

    def __test_mutect2_tumor_normal_paired(self):
        shutil.copy(
            f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa'  # ref_fa should be in workdir in runtime
        )
        actual = Mutect2TumorNormalPaired(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.indir}/tumor_sorted.bam',
            normal_bam=f'{self.indir}/normal_sorted.bam'
        )
        expected = f'{self.outdir}/raw.vcf'
        self.assertFileExists(expected, actual)

    def __test_mutect2_tumor_only(self):
        shutil.copy(
            f'{self.indir}/chr9.fa', f'{self.workdir}/chr9.fa'  # ref_fa should be in workdir in runtime
        )
        actual = Mutect2TumorOnly(self.settings).main(
            ref_fa=f'{self.workdir}/chr9.fa',
            tumor_bam=f'{self.indir}/tumor_sorted.bam',
        )
        expected = f'{self.outdir}/raw.vcf'
        self.assertFileExists(expected, actual)

    def __test_snpeff(self):
        actual = SnpEff(self.settings).main(
            vcf=f'{self.indir}/raw.vcf'
        )
        expected = f'{self.outdir}/annotated.vcf'
        self.assertFileExists(expected, actual)
