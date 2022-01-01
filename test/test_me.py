from .setup import TestCase
from gatk_pipeline.annotation import SnpEff
from gatk_pipeline.trimming import TrimGalore
from gatk_pipeline.variant_calling import Mutect2
from gatk_pipeline.mapping import BwaIndex, BwaMem


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

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

    def __test_mapping(self):
        index = BwaIndex(self.settings).main(
            fna=f'{self.indir}/chr9.fa'
        )
        for sample_name in ['normal', 'tumor']:
            sorted_bam = BwaMem(self.settings).main(
                index=index,
                fq1=f'{self.indir}/{sample_name}.1.fq.gz',
                fq2=f'{self.indir}/{sample_name}.2.fq.gz',
                sample_name=sample_name
            )

    def __test_mutect2(self):
        Mutect2(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa',
            tumor_bam=f'{self.indir}/tumor_sorted.bam',
            normal_bam=f'{self.indir}/normal_sorted.bam'
        )

    def __test_snpeff(self):
        annotated_vcf = SnpEff(self.settings).main(
            vcf=f'{self.indir}/raw.vcf'
        )
        self.assertFileExists(f'{self.outdir}/annotated.vcf', annotated_vcf)
