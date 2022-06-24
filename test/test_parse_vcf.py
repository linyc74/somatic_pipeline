from .setup import TestCase
import numpy as np
import pandas as pd
from somatic_pipeline.parse_vcf import ParseVcf, GetInfoIDToDescription, \
    VcfLineToRow, UnrollSnpEffAnnotation, UnrollVEPAnnotation, GetAllColumns


class TestParseVcf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_snpeff(self):
        ParseVcf(self.settings).main(
            vcf=f'{self.indir}/snpeff.vcf'
        )
        self.assertDataFrameEqual(
            first=pd.read_csv(f'{self.indir}/snpeff.csv'),
            second=pd.read_csv(f'{self.outdir}/variants.csv'),
        )

    def test_vep(self):
        ParseVcf(self.settings).main(
            vcf=f'{self.indir}/vep.vcf'
        )
        self.assertDataFrameEqual(
            first=pd.read_csv(f'{self.indir}/vep.csv'),
            second=pd.read_csv(f'{self.outdir}/variants.csv'),
        )


class TestGetInfoIDToDescription(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        vcf_header = '''\
##SKIP
##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests.">
##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Number of reads with unique start and mate end positions for each alt at a variant site">
##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline variants">
##INFO=<ID=MBQ,Number=R,Type=Integer,Description="median base quality by allele">
##SKIP
##SKIP
##SKIP
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name' ">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	normal	tumor'''

        actual = GetInfoIDToDescription(self.settings).main(vcf_header=vcf_header)

        expected = {
            'AS_SB_TABLE': 'Allele-specific forward/reverse read counts for strand bias tests.',
            'AS_UNIQ_ALT_READ_COUNT': 'Number of reads with unique start and mate end positions for each alt at a variant site',
            'GERMQ': 'Phred-scaled quality that alt alleles are not germline variants',
            'MBQ': 'median base quality by allele',
            'ANN': "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name' ",
        }

        self.assertDictEqual(expected, actual)


class TestGetAllColumns(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = GetAllColumns(self.settings).main(
            info_id_to_description={
                'ID1': 'Description 1',
                'ANN': "Functional annotations: 'Allele | Annotation | Gene_Name' ",
                'CSQ': 'Consequence annotations from Ensembl VEP. Format: ALLELE|CONSEQUENCE|IMPACT',
                'ID2': 'Description 2',
            }
        )
        expected = [
            'Chromosome',
            'Position',
            'ID',
            'Ref Allele',
            'Alt Allele',
            'Quality',
            'Filter',

            'Description 1',
            'Description 2',

            'Allele',
            'Annotation',
            'Gene_Name',

            'ALLELE',
            'CONSEQUENCE',
            'IMPACT',
        ]
        self.assertListEqual(expected, actual)


class TestVcfLineToRow(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        vcf_line = 'chr9\t12303\t.\tA\tT\t.\t.\tMBQ=20,20;ANN=T|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|transcript|ENST00000442898.5|protein_coding||c.*2504T>A|||||2172|,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L5|ENSG00000236875|transcript|ENST00000421620.2|unprocessed_pseudogene|2/6|n.70A>T||||||;dbNSFP_CADD_phred=16.54'
        info_id_to_description = {
            'MBQ': 'median base quality by allele',
            'ANN': 'Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\' ',
            'dbNSFP_CADD_phred': 'dbNSFP CADD Phred Score',
        }
        actual = VcfLineToRow(self.settings).main(
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
            'dbNSFP CADD Phred Score': '16.54',

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


class TestUnrollVEPAnnotation(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_mock(self):
        d = {
            "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|BAM_EDIT":
                'A|consequence(1)|impact(1)|symbol(1)|bam_edit(1)|consequence(2)|impact(2)|symbol(2)|bam_edit(2)'
        }
        actual = UnrollVEPAnnotation(self.settings).main(d=d)
        expected = {
            'Allele': 'A',
            'Consequence': 'consequence(1)',
            'IMPACT': 'impact(1)',
            'SYMBOL': 'symbol(1)',
            'BAM_EDIT': 'bam_edit(1)',
        }
        self.assertDictEqual(expected, actual)

    def test_real(self):
        d = {
            "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT":
                'A|intron_variant&non_coding_transcript_variant|MODIFIER|DDX11L5|ENSG00000236875|Transcript|ENST00000421620|unprocessed_pseudogene||2/5||||||||||1||HGNC|HGNC:37106||Ensembl||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|Transcript|ENST00000442898|protein_coding|||||||||||2077|-1||HGNC|HGNC:24361||Ensembl||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|Transcript|ENST00000696149|protein_coding|||||||||||2115|-1||HGNC|HGNC:24361||Ensembl||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|Transcript|ENST00000696150|processed_transcript|||||||||||3511|-1||HGNC|HGNC:24361||Ensembl||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|ENSG00000181404|Transcript|ENST00000696151|nonsense_mediated_decay|||||||||||4163|-1||HGNC|HGNC:24361||Ensembl||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|NM_001378090.1|protein_coding|||||||||||2115|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|NM_182905.6|protein_coding|||||||||||2115|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|intron_variant&non_coding_transcript_variant|MODIFIER|DDX11L5|100287596|Transcript|NR_051986.1|transcribed_pseudogene||1/2||||||||||1||EntrezGene|HGNC:37106||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517659.2|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517660.2|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517662.3|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517663.3|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517664.3|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517665.2|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_011517666.2|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_017014169.1|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_017014170.1|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_017014171.1|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_017014172.1|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_017014173.1|protein_coding|||||||||||4332|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_024447369.1|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XM_024447370.1|protein_coding|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XR_002956737.1|misc_RNA|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XR_002956738.1|misc_RNA|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XR_002956739.1|misc_RNA|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XR_002956740.1|misc_RNA|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XR_002956741.1|misc_RNA|||||||||||2683|-1||EntrezGene|HGNC:24361||RefSeq||G|G|,A|downstream_gene_variant|MODIFIER|WASHC1|100287171|Transcript|XR_929142.2|misc_RNA|||||||||||2077|-1||EntrezGene|HGNC:24361||RefSeq||G|G|'
        }
        actual = UnrollVEPAnnotation(self.settings).main(d=d)
        expected = {
            'Allele': 'A',
            'Consequence': 'intron_variant&non_coding_transcript_variant',
            'IMPACT': 'MODIFIER',
            'SYMBOL': 'DDX11L5',
            'Gene': 'ENSG00000236875',
            'Feature_type': 'Transcript',
            'Feature': 'ENST00000421620',
            'BIOTYPE': 'unprocessed_pseudogene',
            'EXON': '',
            'INTRON': '2/5',
            'HGVSc': '',
            'HGVSp': '',
            'cDNA_position': '',
            'CDS_position': '',
            'Protein_position': '',
            'Amino_acids': '',
            'Codons': '',
            'Existing_variation': '',
            'DISTANCE': '',
            'STRAND': '1',
            'FLAGS': '',
            'SYMBOL_SOURCE': 'HGNC',
            'HGNC_ID': 'HGNC:37106',
            'REFSEQ_MATCH': '',
            'SOURCE': 'Ensembl',
            'REFSEQ_OFFSET': '',
            'GIVEN_REF': 'G',
            'USED_REF': 'G',
            'BAM_EDIT': ',A',
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
