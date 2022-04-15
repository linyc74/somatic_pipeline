import pandas as pd
from typing import Dict, Any, List
from .template import Processor, Settings


class ParseSnpEffVcf(Processor):

    LOG_INTERVAL = 10000  # variants

    vcf: str

    vcf_header: str
    info_id_to_description: Dict[str, str]
    columns: List[str]
    data: List[Dict[str, Any]]  # each dict is a row (i.e. variant)

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.vcf_line_to_row = SnpEffVcfLineToRow(self.settings).main

    def main(self, vcf: str):
        self.vcf = vcf

        self.logger.info(msg='Start parsing annotated VCF')
        self.set_vcf_header()
        self.set_info_id_to_description()
        self.set_columns()
        self.process_vcf_data()
        self.save_csv()

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

    def set_columns(self):
        info_descriptions = list(self.info_id_to_description.values())
        self.columns = [
            'Chromosome',
            'Position',
            'ID',
            'Ref Allele',
            'Alt Allele',
            'Quality',
            'Filter',
        ] + info_descriptions

    def process_vcf_data(self):
        n = 0
        self.data = []
        with open(self.vcf) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                n += 1
                if n % self.LOG_INTERVAL == 0:
                    self.logger.debug(msg=f'{n} variants parsed')

                row = self.vcf_line_to_row(
                    vcf_line=line,
                    info_id_to_description=self.info_id_to_description)
                self.data.append(row)

    def save_csv(self):
        df = pd.DataFrame(self.data, columns=self.columns)
        df.to_csv(f'{self.outdir}/variants.csv', index=False)


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


class SnpEffVcfLineToRow(Processor):

    vcf_line: str
    info_id_to_description: Dict[str, str]

    vcf_info: str
    row: Dict[str, Any]

    def __init__(self, settings: Settings):
        super().__init__(settings)
        self.unroll_snpeff_annotation = UnrollSnpEffAnnotation(self.settings).main

    def main(
            self,
            vcf_line: str,
            info_id_to_description: Dict[str, str]) -> Dict[str, Any]:

        self.vcf_line = vcf_line
        self.info_id_to_description = info_id_to_description

        self.unpack_line_and_set_vcf_info()
        self.parse_vcf_info()
        self.row = self.unroll_snpeff_annotation(self.row)

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
            if '=' not in item:
                continue

            id_, val = item.split('=')
            description = self.info_id_to_description.get(id_, None)
            if description is not None:
                self.row[description] = val


class UnrollSnpEffAnnotation(Processor):

    LEFT_STRIP = "Functional annotations: '"
    RIGHT_STRIP = "' "

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
        keys = key[len(self.LEFT_STRIP):-len(self.RIGHT_STRIP)].split(' | ')
        vals = val.split('|')
        new_dict = {
            k: v for k, v in zip(keys, vals)
        }
        self.d.update(new_dict)
