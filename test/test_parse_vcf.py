from .setup import TestCase
import numpy as np
import pandas as pd
from somatic_pipeline.parse_vcf import ParseMutect2SnpEffVcf, GetInfoIDToDescription, \
    Mutect2SnpEffVcfLineToRow, UnrollSnpEffAnnotation


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


class TestGetInfoIDToDescription(TestCase):

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
##MutectVersion=2.2
##SnpEffVersion="5.0e (build 2021-03-09 06:01), by Pablo Cingolani"
##SnpEffCmd="SnpEff  -htmlStats test/outdir/snpEff_summary.html GRCh38.99 test/test_me/raw.vcf "
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	normal	tumor'''
        expected = {
            'AS_SB_TABLE': 'Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.',
            'AS_UNIQ_ALT_READ_COUNT': 'Number of reads with unique start and mate end positions for each alt at a variant site',
            'GERMQ': 'Phred-scaled quality that alt alleles are not germline variants',
            'MBQ': 'median base quality by allele',
            'ANN': 'Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\' ',
        }
        actual = GetInfoIDToDescription(self.settings).main(
            vcf_header=vcf_header
        )
        self.assertDictEqual(expected, actual)


class TestMutect2SnpEffVcfLineToRow(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        vcf_line = 'chr9\t12303\t.\tA\tT\t.\t.\tMBQ=20,20;ANN=T|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|transcript|ENST00000442898.5|protein_coding||c.*2504T>A|||||2172|,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L5|ENSG00000236875|transcript|ENST00000421620.2|unprocessed_pseudogene|2/6|n.70A>T||||||'
        info_id_to_description = {
            'MBQ': 'median base quality by allele',
            'ANN': 'Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\' ',
        }
        actual = Mutect2SnpEffVcfLineToRow(self.settings).main(
            vcf_line=vcf_line,
            info_id_to_description=info_id_to_description
        )
        expected = {
            'Chromosome': 'chr9',
            'Position': '12303',
            'ID': '.',
            'Ref Allele': 'A',
            'Alt Allele': 'T',
            'Quality': '.',
            'Filter': '.',
            'median base quality by allele': '20,20',
            'Allele': 'T',
            'Annotation': 'downstream_gene_variant',
            'Annotation_Impact': 'MODIFIER',
            'Gene_Name': 'WASHC1',
            'Gene_ID': 'ENSG00000181404',
            'Feature_Type': 'transcript',
            'Feature_ID': 'ENST00000442898.5',
            'Transcript_BioType': 'protein_coding',
            'Rank': '',
            'HGVS.c': 'c.*2504T>A',
            'HGVS.p': '',
            'cDNA.pos / cDNA.length': '',
            'CDS.pos / CDS.length': '',
            'AA.pos / AA.length': '',
            'Distance': '2172',
            'ERRORS / WARNINGS / INFO': ',T'
        }
        self.assertDictEqual(expected, actual)


class TestUnrollSnpEffAnnotation(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        d = {
            "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ":
                'T|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|transcript|ENST00000442898.5|protein_coding||c.*2504T>A|||||2172|,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L5|ENSG00000236875|transcript|ENST00000421620.2|unprocessed_pseudogene|2/6|n.70A>T||||||',
        }
        actual = UnrollSnpEffAnnotation(self.settings).main(d=d)
        expected = {
            'Allele': 'T',
            'Annotation': 'downstream_gene_variant',
            'Annotation_Impact': 'MODIFIER',
            'Gene_Name': 'WASHC1',
            'Gene_ID': 'ENSG00000181404',
            'Feature_Type': 'transcript',
            'Feature_ID': 'ENST00000442898.5',
            'Transcript_BioType': 'protein_coding',
            'Rank': '',
            'HGVS.c': 'c.*2504T>A',
            'HGVS.p': '',
            'cDNA.pos / cDNA.length': '',
            'CDS.pos / CDS.length': '',
            'AA.pos / AA.length': '',
            'Distance': '2172',
            'ERRORS / WARNINGS / INFO': ',T',
        }
        self.assertDictEqual(expected, actual)


class TestPandasDataStructure(TestCase):
    """
    To understand how to append rows without using df.append(), which is extremely slow
    """

    def test_list_of_dict(self):
        df1 = pd.DataFrame([
            {'A': 1, 'B': 1},  # row 1
            {'A': 2, 'C': 2},  # row 2
            {'D': 3},          # row 3
            {'A': 4},          # row 4
        ])

        df2 = pd.DataFrame(columns=[
            'A', 'B', 'C', 'D'
        ], data=[
            [1, 1, np.nan, np.nan],
            [2, np.nan, 2, np.nan],
            [np.nan, np.nan, np.nan, 3],
            [4, np.nan, np.nan, np.nan]
        ])

        self.assertDataFrameEqual(df1, df2)
