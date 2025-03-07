import pandas as pd
from os.path import exists, dirname
from typing import Dict, Any, List, Tuple, Optional
from .tools import edit_fpath
from .template import Processor, Settings


class Vcf2Csv(Processor):

    LOG_INTERVAL = 10000  # variants

    vcf: str
    dstdir: Optional[str]

    vcf_header: str
    info_id_to_description: Dict[str, str]
    all_columns: List[str]
    output_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.vcf_line_to_row = VcfLineToRow(self.settings).main
        self.save_data_to_csv = SaveDataToCsv(self.settings).main

    def main(self, vcf: str, dstdir: Optional[str] = None) -> str:
        self.vcf = vcf
        self.dstdir = dstdir

        self.logger.info(msg=f'Start parsing VCF: {self.vcf}')
        self.set_vcf_header()
        self.set_info_id_to_description()
        self.set_all_columns()
        self.set_output_csv()
        self.process_vcf_data()

        return self.output_csv

    def set_vcf_header(self):
        self.vcf_header = ''
        with open(self.vcf) as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                self.vcf_header += line

    def set_info_id_to_description(self):
        self.info_id_to_description = GetInfoIDToDescription(self.settings).main(
            vcf_header=self.vcf_header)

    def set_all_columns(self):
        self.all_columns = GetAllColumns(self.settings).main(
            info_id_to_description=self.info_id_to_description)

    def set_output_csv(self):
        self.output_csv = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='.csv',
            dstdir=dirname(self.vcf) if self.dstdir is None else self.dstdir
        )

    def process_vcf_data(self):
        n = 0
        data: List[Dict[str, Any]]  # each dict is a row (i.e. variant)
        data = []
        with open(self.vcf) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                row = self.__line_to_row(line)
                data.append(row)

                n += 1
                if n % self.LOG_INTERVAL == 0:
                    self.logger.debug(msg=f'{n} variants parsed')
                    self.__to_csv(data=data)
                    data = []  # clear up data
            self.__to_csv(data=data)  # last partial chunk of data

    def __line_to_row(self, line: str) -> Dict[str, Any]:
        return self.vcf_line_to_row(
            vcf_line=line,
            info_id_to_description=self.info_id_to_description)

    def __to_csv(self, data: List[Dict[str, Any]]):
        self.save_data_to_csv(
            data=data,
            all_columns=self.all_columns,
            csv=self.output_csv
        )


class GetInfoIDToDescription(Processor):

    vcf_header: str

    info_lines: List[str]
    id_to_description: Dict[str, str]

    def main(self, vcf_header: str) -> Dict[str, str]:
        self.vcf_header = vcf_header

        self.set_info_lines()
        self.set_id_to_description()

        return self.id_to_description

    def set_info_lines(self):
        self.info_lines = []
        for line in self.vcf_header.splitlines():
            if line.startswith('##INFO'):
                self.info_lines.append(line)

    def set_id_to_description(self):
        self.id_to_description = {}
        for line in self.info_lines:
            self.process_one(info_line=line)

    def process_one(self, info_line: str):
        """
        ##INFO=<ID=MBQ,Number=R,Type=Integer,Description="median base quality by allele">

        id_ = 'MBQ'
        description = 'median base quality by allele'
        """
        id_ = info_line.split('INFO=<ID=')[1].split(',')[0]
        description = info_line.split(',Description="')[1].split('">')[0]
        self.id_to_description[id_] = description


class GetAllColumns(Processor):

    BASE_COLUMNS = [
        'Chromosome',
        'Position',
        'ID',
        'Ref Allele',
        'Alt Allele',
        'Quality',
        'Filter',
    ]

    info_id_to_description: Dict[str, str]

    columns: List[str]

    def main(
            self,
            info_id_to_description: Dict[str, str]) -> List[str]:

        self.info_id_to_description = info_id_to_description

        self.init_columns()
        self.add_info_descriptions()
        self.unroll_vep_columns()

        return self.columns

    def init_columns(self):
        self.columns = self.BASE_COLUMNS.copy()

    def add_info_descriptions(self):
        for description in self.info_id_to_description.values():
            self.columns.append(description)

    def unroll_vep_columns(self):
        """
        "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT"

        ['Allele', 'Consequence', 'IMPACT']
        """
        left_strip = 'Consequence annotations from Ensembl VEP. Format: '
        sep = '|'
        for i, column in enumerate(self.columns):
            if column.startswith(left_strip):
                self.columns.pop(i)
                new_columns = column.lstrip(left_strip).split(sep)
                self.columns += new_columns
                break


class VcfLineToRow(Processor):

    vcf_line: str
    info_id_to_description: Dict[str, str]

    vcf_info: str
    row: Dict[str, Any]

    def main(
            self,
            vcf_line: str,
            info_id_to_description: Dict[str, str]) -> Dict[str, Any]:

        self.vcf_line = vcf_line
        self.info_id_to_description = info_id_to_description

        self.unpack_line_and_set_vcf_info()
        self.parse_vcf_info()
        self.row = UnrollVEPAnnotation(self.settings).main(row=self.row)

        return self.row

    def unpack_line_and_set_vcf_info(self):
        chromosome, position, id_, ref_allele, alt_allele, quality, filter_, info = \
            self.vcf_line.strip().split('\t')[:8]

        self.row = {
            'Chromosome': chromosome,
            'Position': position,
            'ID': id_,
            'Ref Allele': ref_allele,
            'Alt Allele': alt_allele,
            'Quality': quality,
            'Filter': filter_,
        }
        self.vcf_info = info

    def parse_vcf_info(self):
        items = self.vcf_info.split(';')
        for item in items:
            if '=' in item:
                key, val = self.__split(item)
                description = self.info_id_to_description.get(key, None)

            else:  # there is no '=', this is a flag, without value
                key, val = item, True
                description = self.info_id_to_description.get(item, None)

            if description is not None:
                self.row[description] = val

    def __split(self, item: str) -> Tuple[str, str]:
        p = item.index('=')
        id_ = item[0:p]
        val = item[p+1:]
        return id_, val


class UnrollVEPAnnotation(Processor):

    VEP_PREFIX = 'Consequence annotations from Ensembl VEP. Format: '
    # https://gdc.cancer.gov/content/most-frequent-mutations-table-vep-impact-score-which-algorithm-vep-gdc-using-determine-%E2%80%9Ch-or-%E2%80%9Cm%E2%80%9D
    IMPACT_ORDER = [
        'HIGH',
        'MODERATE',
        'LOW',
        'MODIFIER',
    ]
    CONSEQUENCE_ORDER = [
        # amino acid change
        'frameshift_variant',
        'stop_gained',
        'missense_variant',
        'inframe_deletion',
        'inframe_insertion',
        'stop_lost',
        'start_lost',
        'stop_retained_variant',
        'start_retained_variant',
        'incomplete_terminal_codon_variant',
        'coding_sequence_variant',
        'protein_altering_variant',

        # silent
        'synonymous_variant',

        # splice
        'splice_polypyrimidine_tract_variant',
        'splice_donor_variant',
        'splice_acceptor_variant',
        'splice_donor_region_variant',
        'splice_donor_5th_base_variant',
        'splice_region_variant',

        # UTR
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',

        # regulatory
        'upstream_gene_variant',
        'downstream_gene_variant',
        'regulatory_region_ablation',
        'regulatory_region_variant',
        'TF_binding_site_variant',
        'TFBS_ablation',

        # intron and intergenic
        'intron_variant',
        'intergenic_variant',

        # miRNA and lincRNA
        'mature_miRNA_variant',
        'NMD_transcript_variant',
        'non_coding_transcript_exon_variant',
        'non_coding_transcript_variant',
    ]

    row: Dict[str, str]

    vep_key: Optional[str]
    annotations: List[str]
    fields: List[str]
    values: List[str]

    def main(self, row: Dict[str, str]) -> Dict[str, str]:
        self.row = row.copy()

        self.set_vep_key()
        if self.vep_key is None:  # no VEP annotation
            return self.row

        self.set_annotations()
        self.annotations = sorted(self.annotations, key=self.get_order)
        self.set_fields_and_values()
        self.pop_and_replace_with_new_dict()

        return self.row

    def set_vep_key(self):
        self.vep_key = None
        for key in self.row.keys():
            if key.startswith(self.VEP_PREFIX):
                self.vep_key = key
                break

    def set_annotations(self):
        full_annotation = self.row[self.vep_key]
        self.annotations = full_annotation.split(',')

    def get_order(self, annotation: str) -> Tuple[int, int, int]:
        consequence = annotation.split('|')[1].split('&')[0]
        impact = annotation.split('|')[2]
        symbol = annotation.split('|')[3]

        s, i, c = 999, 999, 999  # last in order by default

        if symbol != '':  # has symbol
            s = 0

        if impact in self.IMPACT_ORDER:
            i = self.IMPACT_ORDER.index(impact)

        if consequence in self.CONSEQUENCE_ORDER:
            c = self.CONSEQUENCE_ORDER.index(consequence)

        return s, i, c

    def set_fields_and_values(self):
        self.fields = self.vep_key[len(self.VEP_PREFIX):].split('|')
        self.values = self.annotations[0].split('|')

    def pop_and_replace_with_new_dict(self):
        self.row.pop(self.vep_key)
        new_dict = {
            k: v for k, v in zip(self.fields, self.values)
        }
        self.row.update(new_dict)


class SaveDataToCsv(Processor):

    data: List[Dict[str, Any]]  # each dict is a row (i.e. variant)
    all_columns: List[str]
    csv: str

    def main(
            self,
            data: List[Dict[str, Any]],
            all_columns: List[str],
            csv: str):

        self.data = data
        self.all_columns = all_columns
        self.csv = csv

        self.write_to_csv()

    def write_to_csv(self):
        if not exists(self.csv):
            header = True
        else:
            header = False

        pd.DataFrame(
            data=self.data,
            columns=self.all_columns
        ).to_csv(
            self.csv,
            mode='a',
            header=header,
            index=False
        )
