import pandas as pd
from os.path import exists
from typing import Dict, Any, List, Tuple
from .template import Processor, Settings


class ParseVcf(Processor):

    LOG_INTERVAL = 10000  # variants

    vcf: str

    vcf_header: str
    info_id_to_description: Dict[str, str]
    all_columns: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.vcf_line_to_row = VcfLineToRow(self.settings).main
        self.save_data_to_csv = SaveDataToCsv(self.settings).main

    def main(self, vcf: str):
        self.vcf = vcf

        self.logger.info(msg='Start parsing annotated VCF')
        self.set_vcf_header()
        self.set_info_id_to_description()
        self.set_all_columns()
        self.process_vcf_data()

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
            csv=f'{self.outdir}/variants.csv'
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
        self.unroll_snpeff_columns()
        self.unroll_vep_columns()

        return self.columns

    def init_columns(self):
        self.columns = self.BASE_COLUMNS.copy()

    def add_info_descriptions(self):
        for description in self.info_id_to_description.values():
            self.columns.append(description)

    def unroll_snpeff_columns(self):
        """
        "Functional annotations: 'Allele | Annotation | Gene_Name' "

        ['Allele', 'Annotation', 'Gene_Name']
        """
        left_strip = "Functional annotations: '"
        right_strip = "' "
        sep = ' | '
        for i, column in enumerate(self.columns):
            if column.startswith(left_strip):
                self.columns.pop(i)
                new_columns = column.lstrip(left_strip).rstrip(right_strip).split(sep)
                self.columns += new_columns
                break

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
        self.unroll_annotation()

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
                id_, val = self.__split(item)
                description = self.info_id_to_description.get(id_, None)
                if description is not None:
                    self.row[description] = val

    def __split(self, item: str) -> Tuple[str, str]:
        p = item.index('=')
        id_ = item[0:p]
        val = item[p+1:]
        return id_, val

    def unroll_annotation(self):
        for unroller in [
            UnrollSnpEffAnnotation(self.settings),
            UnrollVEPAnnotation(self.settings)
        ]:
            self.row = unroller.main(self.row)


class UnrollAnnotation(Processor):

    LEFT_STRIP: str
    RIGHT_STRIP: str
    KEY_SEP: str
    VALUE_SEP: str

    d: Dict[str, str]

    def main(self, d: Dict[str, str]) -> Dict[str, str]:
        self.d = d.copy()

        keys = list(self.d.keys())
        for key in keys:
            if key.startswith(self.LEFT_STRIP):
                val = self.d.pop(key)
                self.unroll(key, val)

        return self.d

    def unroll(self, key: str, val: str):
        keys = key.lstrip(self.LEFT_STRIP).rstrip(self.RIGHT_STRIP).split(self.KEY_SEP)
        vals = val.split(self.VALUE_SEP)
        new_dict = {
            k: v for k, v in zip(keys, vals)
        }
        self.d.update(new_dict)


class UnrollSnpEffAnnotation(UnrollAnnotation):

    LEFT_STRIP = "Functional annotations: '"
    RIGHT_STRIP = "' "
    KEY_SEP = ' | '
    VALUE_SEP = '|'


class UnrollVEPAnnotation(UnrollAnnotation):

    LEFT_STRIP = 'Consequence annotations from Ensembl VEP. Format: '
    RIGHT_STRIP = ''
    KEY_SEP = '|'
    VALUE_SEP = '|'


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
