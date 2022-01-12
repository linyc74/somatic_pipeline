from .setup import TestCase
import pandas as pd
from gatk_pipeline.parse_vcf import ParseMutect2SnpEffVcf, GetSnpEffAnnotationKeys, GetMutect2InfoKeyToName


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


class TestGetMutect2InfoKeyToName(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        vcf_header = '''\
##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2 --normal-sample normal --output test/outdir/raw.vcf ...
##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.">
##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Number of reads with unique start and mate end positions for each alt at a variant site">
##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline variants">
##INFO=<ID=MBQ,Number=R,Type=Integer,Description="median base quality by allele">
##MutectVersion=2.2'''

        expected = {
            'AS_SB_TABLE': 'Mutect2 Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.',
            'AS_UNIQ_ALT_READ_COUNT': 'Mutect2 Number of reads with unique start and mate end positions for each alt at a variant site',
            'GERMQ': 'Mutect2 Phred-scaled quality that alt alleles are not germline variants',
            'MBQ': 'Mutect2 median base quality by allele',

        }
        actual = GetMutect2InfoKeyToName(self.settings).main(
            vcf_header=vcf_header
        )
        self.assertDictEqual(expected, actual)


class TestGetSnpEffAnnotationKeys(TestCase):

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

        actual = GetSnpEffAnnotationKeys(self.settings).main(
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
