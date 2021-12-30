from .setup import TestCase
import pandas as pd
from gatk_pipeline.parse_vcf import ParseMutect2SnpEffVcf, GetSnpEffAnnotationKeysFromVcfHeader


class TestParseMutect2SnpEffVcf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        ParseMutect2SnpEffVcf(self.settings).main(
            vcf=f'{self.indir}/annotated.vcf'
        )
        self.assertDataFrameEqual(
            first=pd.read_csv(f'{self.indir}/variants.csv'),
            second=pd.read_csv(f'{self.outdir}/variants.csv'),
        )


class TestGetSnpEffAnnotationKeysFromVcfHeader(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        vcf_header = '''\
##SnpEffVersion="5.0e (build 2021-03-09 06:01), by Pablo Cingolani"
##SnpEffCmd="SnpEff  -htmlStats test/outdir/snpEff_summary.html GRCh38.99 test/test_me/raw.vcf "
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	normal	tumor'''

        actual = GetSnpEffAnnotationKeysFromVcfHeader(self.settings).main(
            vcf_header=vcf_header
        )
        expected = [
            'SnpEff Allele',
            'SnpEff Annotation',
            'SnpEff Annotation_Impact',
            'SnpEff Gene_Name',
            'SnpEff Gene_ID',
            'SnpEff Feature_Type',
            'SnpEff Feature_ID',
            'SnpEff Transcript_BioType',
            'SnpEff Rank',
            'SnpEff HGVS.c',
            'SnpEff HGVS.p',
            'SnpEff cDNA.pos / cDNA.length',
            'SnpEff CDS.pos / CDS.length',
            'SnpEff AA.pos / AA.length',
            'SnpEff Distance',
            'SnpEff ERRORS / WARNINGS / INFO',
        ]
        self.assertListEqual(expected, actual)
