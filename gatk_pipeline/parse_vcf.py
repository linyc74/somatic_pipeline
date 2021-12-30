import pandas as pd
from typing import Dict, Any, List
from .template import Processor, Settings


class ParseMutect2SnpEffVcf(Processor):

    vcf: str

    vcf_header: str
    snpeff_annotation_keys: List[str]
    df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, vcf: str):
        self.vcf = vcf

        self.set_vcf_header()
        self.set_snpeff_annotation_keys()
        self.process_vcf_data()
        self.save_csv()

    def set_vcf_header(self):
        self.vcf_header = ''
        with open(self.vcf) as fh:
            for line in fh:
                if line.startswith('##'):
                    self.vcf_header += line

    def set_snpeff_annotation_keys(self):
        self.snpeff_annotation_keys = GetSnpEffAnnotationKeysFromVcfHeader(self.settings).main(
            vcf_header=self.vcf_header)

    def process_vcf_data(self):
        self.df = pd.DataFrame()
        line_to_row = Mutect2SnpEffVcfLineToRow(self.settings).main

        with open(self.vcf) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                row = line_to_row(
                    vcf_line=line,
                    snpeff_annotation_keys=self.snpeff_annotation_keys)

                self.df = self.df.append(row, ignore_index=True)

    def save_csv(self):
        self.df.to_csv(f'{self.outdir}/variants.csv', index=False)


class GetSnpEffAnnotationKeysFromVcfHeader(Processor):

    KEY_PREFIX = 'SnpEff '

    vcf_header: str
    keys: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, vcf_header: str):
        self.vcf_header = vcf_header

        self.keys = []

        for line in self.vcf_header.splitlines():
            if line.startswith('##INFO=<ID=ANN'):
                self.process_(annotation_line=line)

        self.add_prefix()

        return self.keys

    def process_(self, annotation_line: str):
        """
        ##INFO=<ID=ANN,...,Description="Functional annotations: 'Allele | Annotation | ... | Distance | ERRORS / WARNINGS / INFO' ">
        """
        after_this = 'Description="Functional annotations:'
        before_this = '">'

        line = annotation_line
        middle = line.split(after_this)[1].split(before_this)[0].strip()
        middle = middle[1:-1]  # Remove ' and ' at the beginning and end

        self.keys += middle.split(' | ')

    def add_prefix(self):
        for i, key in enumerate(self.keys):
            self.keys[i] = self.KEY_PREFIX + key


class Mutect2SnpEffVcfLineToRow(Processor):

    vcf_line: str
    snpeff_annotation_keys: List[str]

    vcf_info: str
    row: Dict[str, Any]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            vcf_line: str,
            snpeff_annotation_keys: List[str]) -> Dict[str, Any]:

        self.vcf_line = vcf_line
        self.snpeff_annotation_keys = snpeff_annotation_keys

        self.unpack_line()
        self.parse_snpeff_annotation_from_vcf_info()

        return self.row

    def unpack_line(self):
        chromosome, position, id_, ref_allele, \
        alt_allele, quality, filter_, info = \
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

    def parse_snpeff_annotation_from_vcf_info(self):
        ann_str = None
        items = self.vcf_info.split(';')
        for item in items:
            if item.startswith('ANN='):
                ann_str = item[len('ANN='):]

        if ann_str is None:
            return

        ann_values = ann_str.split('|')

        ann_dict = {
            k: v for k, v in
            zip(self.snpeff_annotation_keys, ann_values)
        }

        self.row.update(ann_dict)
