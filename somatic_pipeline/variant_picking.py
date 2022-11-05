import gzip
from os.path import basename
from typing import List, IO, Optional, Dict, Any, Tuple
from .template import Processor
from .constant import TUMOR, NORMAL


class VcfParser:

    header: str
    columns: List[str]

    __fh: IO

    def __init__(self, vcf: str):
        if vcf.endswith('.gz'):
            self.__fh = gzip.open(vcf, 'rt')  # rt: read text
        else:
            self.__fh = open(vcf, 'r')
        self.__set_header()

    def __set_header(self):
        header_lines = []
        for line in self.__fh:
            header_lines.append(line.strip())

            if line.startswith('##'):
                continue
            else:
                assert line.startswith('#')
                break

        self.header = '\n'.join(header_lines)
        self.columns = header_lines[-1][1:].split('\t')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        r = self.next()
        if r is not None:
            return r
        else:
            raise StopIteration

    def next(self) -> Optional[Dict[str, str]]:
        line = self.__fh.readline()

        if line == '':  # end of file
            return None

        values = line.strip().split('\t')

        assert len(self.columns) == len(values)

        return {c: v for c, v in zip(self.columns, values)}

    def close(self):
        self.__fh.close()


class VcfWriter:

    header: Optional[str]
    columns: List[str]

    __fh: IO

    def __init__(self, vcf: str):
        self.__fh = open(vcf, 'w')
        self.header = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def write_header(self, header: str):
        assert self.header is None
        self.header = header.strip()
        self.__assert_header_format()
        self.__set_columns()
        self.__fh.write(self.header + '\n')

    def __assert_header_format(self):
        header_lines = self.header.splitlines()

        for line in header_lines[:-1]:
            assert line.startswith('##')

        last_line = header_lines[-1]
        assert last_line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

    def __set_columns(self):
        last_line = self.header.splitlines()[-1]
        self.columns = last_line[1:].split('\t')

    def write(self, variant: Dict[str, Any]):
        assert self.header is not None  # header must have been written
        assert set(variant.keys()) == set(self.columns)

        # values need to follow the order of self.columns
        values = [variant[c] for c in self.columns]

        self.__fh.write('\t'.join(map(str, values)) + '\n')

    def close(self):
        self.__fh.close()


class VariantPicking(Processor):

    VARIANT_KEY_COLUMNS = ['CHROM', 'POS', 'REF', 'ALT']
    INFO_CALL_ID = 'CALL'
    INFO_CALL_DESCRIPTION = 'Variant callers detecting this variant'

    ref_fa: str
    vcfs: List[str]
    min_num_callers: int

    variant_to_callers: Dict[Tuple, List[str]]
    header: str
    out_vcf: str

    def main(
            self,
            ref_fa: str,
            vcfs: List[str],
            min_num_callers: int) -> str:

        self.ref_fa = ref_fa
        self.vcfs = vcfs
        self.min_num_callers = min_num_callers

        self.collect_variant_dict()
        self.set_header()
        self.write_vcf()

        return self.out_vcf

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

    def set_header(self):
        contig_lines = BuildHeaderContigLines(self.settings).main(self.ref_fa)
        self.header = f'''\
##fileformat=VCFv4.2
{contig_lines}
##INFO=<ID={self.INFO_CALL_ID},Number=.,Type=String,Description="{self.INFO_CALL_DESCRIPTION}">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{NORMAL}\t{TUMOR}'''

    def write_vcf(self):
        self.out_vcf = f'{self.outdir}/variants.vcf'
        with VcfWriter(self.out_vcf) as writer:
            writer.write_header(self.header)
            for variant, callers in self.variant_to_callers.items():
                if len(callers) >= self.min_num_callers:
                    v = {c: '.' for c in writer.columns}  # empty dict
                    v['CHROM'] = variant[0]
                    v['POS'] = variant[1]
                    v['REF'] = variant[2]
                    v['ALT'] = variant[3]
                    c = ','.join(callers)
                    v['INFO'] = f'CALL={c}'
                    writer.write(v)


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
