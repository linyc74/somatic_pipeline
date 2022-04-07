from .setup import TestCase
from somatic_pipeline.somatic_pipeline import SomaticPipeline


class TestSomaticPipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_tumor_normal_paired(self):
        SomaticPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa.gz',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=f'{self.indir}/normal.1.fq.gz',
            normal_fq2=f'{self.indir}/normal.2.fq.gz',
            read_aligner='bwa',
            variant_caller='mutect2',
            exome_target_bed=f'{self.indir}/chr9-exome-probes.bed',
            cnvkit_annotate_txt=f'{self.indir}/chr9-refFlat.txt',
            panel_of_normal_vcf=None,
            bqsr_known_variant_vcf=f'{self.indir}/known-variants.vcf',
            discard_bam=True,
            skip_mark_duplicates=False,
            skip_variant_calling=False,
            skip_cnv=False
        )

    def test_tumor_only(self):
        SomaticPipeline(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa.gz',
            tumor_fq1=f'{self.indir}/tumor.1.fq.gz',
            tumor_fq2=f'{self.indir}/tumor.2.fq.gz',
            normal_fq1=None,
            normal_fq2=None,
            read_aligner='bwa',
            variant_caller='mutect2',
            exome_target_bed=None,
            cnvkit_annotate_txt=None,
            panel_of_normal_vcf=f'{self.indir}/22_0406_twb_snp_pon.vcf',
            bqsr_known_variant_vcf=f'{self.indir}/known-variants.vcf',
            discard_bam=True,
            skip_mark_duplicates=False,
            skip_variant_calling=False,
            skip_cnv=False
        )
