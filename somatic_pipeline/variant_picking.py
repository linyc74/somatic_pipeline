import pandas as pd
from os.path import basename
from typing import List, Dict, Tuple
from .tools import VcfParser
from .template import Processor


class VariantPicking(Processor):

    VARIANT_KEY_COLUMNS = ['CHROM', 'POS', 'REF', 'ALT']
    VCF_COLUMNS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    INFO_CALL_ID = 'CALL'
    INFO_CALL_DESCRIPTION = 'Variant callers detecting this variant'

    ref_fa: str
    vcfs: List[str]
    min_snv_callers: int
    min_indel_callers: int

    vcf_header: str
    variant_to_callers: Dict[Tuple, List[str]]
    variant_df: pd.DataFrame
    out_vcf: str

    def main(
            self,
            ref_fa: str,
            vcfs: List[str],
            min_snv_callers: int,
            min_indel_callers: int) -> str:

        self.ref_fa = ref_fa
        self.vcfs = vcfs
        self.min_snv_callers = min_snv_callers
        self.min_indel_callers = min_indel_callers

        self.build_vcf_header()
        self.collect_variant_dict()
        self.build_variant_df()
        self.sort_variant_df()
        self.write_vcf()

        return self.out_vcf

    def build_vcf_header(self):
        contig_lines = BuildHeaderContigLines(self.settings).main(self.ref_fa)
        columns = '\t'.join(self.VCF_COLUMNS)
        self.vcf_header = f'''\
##fileformat=VCFv4.2
{contig_lines}
##INFO=<ID={self.INFO_CALL_ID},Number=.,Type=String,Description="{self.INFO_CALL_DESCRIPTION}">
#{columns}'''

    def collect_variant_dict(self):
        self.variant_to_callers = {}
        for vcf in self.vcfs:
            caller = self.__get_filename(path=vcf)
            with VcfParser(vcf) as parser:
                for variant in parser:
                    tup = tuple(variant[k] for k in self.VARIANT_KEY_COLUMNS)
                    self.variant_to_callers.setdefault(tup, [])
                    self.variant_to_callers[tup].append(caller)

    def __get_filename(self, path: str) -> str:
        return basename(path).split('.')[0]

    def build_variant_df(self):
        data = []
        for variant, callers in self.variant_to_callers.items():
            chrom, pos, ref, alt = variant

            first_alt = alt.split(',')[0].split('/')[0]
            is_snv = len(ref) == len(first_alt)
            satisfy_snv = len(callers) >= self.min_snv_callers

            is_indel = not is_snv
            satisfy_indel = len(callers) >= self.min_indel_callers

            if (is_snv and satisfy_snv) or (is_indel and satisfy_indel):
                v = {c: '.' for c in self.VCF_COLUMNS}  # empty dict
                v['CHROM'] = chrom
                v['POS'] = int(pos)
                v['REF'] = ref
                v['ALT'] = alt
                c = ','.join(callers)
                v['INFO'] = f'CALL={c}'
                data.append(v)

        self.variant_df = pd.DataFrame(data=data, columns=self.VCF_COLUMNS)

    def sort_variant_df(self):
        chrom_to_order = GetChromToOrder(self.settings).main(self.ref_fa)

        self.variant_df['order'] = self.variant_df['CHROM'].map(chrom_to_order)

        self.variant_df = self.variant_df.sort_values(
            by=['order', 'POS'],
            ascending=[True, True]
        ).drop(
            columns='order'
        )

    def write_vcf(self):
        self.out_vcf = f'{self.workdir}/picked-variants.vcf'

        with open(self.out_vcf, 'w') as writer:
            writer.write(self.vcf_header + '\n')

        self.variant_df.to_csv(
            self.out_vcf,
            mode='a',
            sep='\t',
            header=False,
            index=False)


class BuildHeaderContigLines(Processor):

    ref_fa: str

    contig_id_length: Dict[str, int]
    contig_lines: str

    def main(self, ref_fa: str) -> str:
        self.ref_fa = ref_fa

        self.set_contig_id_to_length()
        self.set_contig_lines()

        return self.contig_lines

    def set_contig_id_to_length(self):
        self.contig_id_length = {}
        with open(self.ref_fa) as fh:
            for line in fh:
                if line.startswith('>'):
                    this_id = line.strip()[1:].split(' ')[0]
                    self.contig_id_length.setdefault(this_id, 0)
                else:
                    self.contig_id_length[this_id] += len(line.strip())

    def set_contig_lines(self):
        lines = [
            f'##contig=<ID={contig_id},length={length}>'
            for contig_id, length in self.contig_id_length.items()
        ]
        self.contig_lines = '\n'.join(lines)


class GetChromToOrder(Processor):

    ref_fa: str

    chrom_to_order: Dict[str, int]

    def main(self, ref_fa: str) -> Dict[str, int]:
        self.ref_fa = ref_fa

        self.chrom_to_order = {}
        current = 0
        with open(self.ref_fa) as fh:
            for line in fh:
                if line.startswith('>'):
                    chrom = line[1:].strip().split(' ')[0]
                    self.chrom_to_order[chrom] = current
                    current += 1

        return self.chrom_to_order
