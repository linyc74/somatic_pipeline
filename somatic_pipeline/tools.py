import os
import gzip
from typing import Optional, List, IO, Dict, Any


def edit_fpath(
        fpath: str,
        old_suffix: str = '',
        new_suffix: str = '',
        dstdir: Optional[str] = None) -> str:

    f = os.path.basename(fpath)
    if old_suffix != '':
        f = f[:-len(old_suffix)]  # strip old suffix
    f += new_suffix

    if dstdir is None:
        dstdir = os.path.dirname(fpath)

    return f'{dstdir}/{f}'


def get_temp_path(
        prefix: str = 'temp',
        suffix: str = '') -> str:

    i = 1
    while True:
        fpath = f'{prefix}{i:03}{suffix}'
        if not os.path.exists(fpath):
            return fpath
        i += 1


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
