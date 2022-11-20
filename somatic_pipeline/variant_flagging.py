from .template import Processor
from .tools import edit_fpath, VcfWriter, VcfParser
from typing import Dict, Any, Optional, Tuple, List


class Criterion:

    def __init__(
            self,
            key: str,
            range: Tuple[float, float],
            equal_max: bool,
            equal_min: bool):

        self.key = key
        self.range = range
        self.equal_max = equal_max
        self.equal_min = equal_min

    def __repr__(self):
        return f"Criterion(key='{self.key}', range={self.range}, equal_max={self.equal_max}, equal_min={self.equal_min})"


class FlagVariants(Processor):

    vcf: str
    variant_flagging_criteria: str

    parser: VcfParser
    writer: VcfWriter
    new_header_lines: List[str]
    flag_to_criterion: Dict[str, Criterion]
    output_vcf: str

    def main(
            self,
            vcf: str,
            variant_flagging_criteria: str) -> str:

        self.vcf = vcf
        self.variant_flagging_criteria = variant_flagging_criteria.replace(' ', '')

        self.open_files()
        self.unpack_variant_flagging_criteria()
        self.write_header()
        self.flag_variants()
        self.close_files()

        return self.output_vcf

    def open_files(self):
        self.parser = VcfParser(self.vcf)
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-flagged.vcf',
            dstdir=self.workdir)
        self.writer = VcfWriter(self.output_vcf)

    def unpack_variant_flagging_criteria(self):
        self.new_header_lines = []
        self.flag_to_criterion = {}

        for item in self.variant_flagging_criteria.split(','):
            flag, criterion = item.split(':')

            self.new_header_lines.append(f'##FILTER=<ID={flag},Description="{criterion}">')

            self.flag_to_criterion[flag] = parse_criterion(s=criterion)

        self.__log()

    def __log(self):
        t = '\n'.join(self.new_header_lines)
        msg = f'Flag variants in "{self.vcf}" with criteria:\n{t}'
        self.logger.info(msg)

    def write_header(self):
        lines = self.parser.header.splitlines()
        lines = lines[0:1] + self.new_header_lines + lines[1:]
        self.writer.write_header('\n'.join(lines))

    def flag_variants(self):
        for variant in self.parser:
            for flag, criterion in self.flag_to_criterion.items():
                variant = flag_variant(
                    variant=variant,
                    flag=flag,
                    criterion=criterion)
            self.writer.write(variant=variant)

    def close_files(self):
        self.parser.close()
        self.writer.close()


def flag_variant(
        variant: Dict[str, Any],
        flag: str,
        criterion: Criterion) -> Dict[str, Any]:

    variant = variant.copy()

    val = get_info_value(variant=variant, key=criterion.key)

    if val is None:
        return variant

    min_, max_ = criterion.range

    if criterion.equal_max:
        less_than_max = val <= max_
    else:
        less_than_max = val < max_

    if criterion.equal_min:
        more_than_min = min_ <= val
    else:
        more_than_min = min_ < val

    if less_than_max and more_than_min:
        variant['FILTER'] += f';{flag}'

    if variant['FILTER'].startswith('.;'):
        variant['FILTER'] = variant['FILTER'][2:]

    return variant


def parse_criterion(s: str) -> Criterion:

    inclusive_min, inclusive_max = False, False

    if '<' in s:

        items = s.split('<')

        if len(items) == 2:
            m, k, M = float('-inf'), items[0], items[1]
        elif len(items) == 3:
            m, k, M = items
        else:
            raise AssertionError

        if k.startswith('='):
            k = k[1:]
            inclusive_min = True

        if M.startswith('='):
            M = M[1:]
            inclusive_max = True

    elif '>' in s:
        items = s.split('>')

        if len(items) == 2:
            M, k, m = float('inf'), items[0], items[1]
        elif len(items) == 3:
            M, k, m = items
        else:
            raise AssertionError

        if k.startswith('='):
            k = k[1:]
            inclusive_max = True

        if m.startswith('='):
            m = m[1:]
            inclusive_min = True

    else:
        raise AssertionError

    min_ = m if m is None else float(m)
    max_ = M if M is None else float(M)

    ret = Criterion(
        key=k,
        range=(min_, max_),
        equal_max=inclusive_max,
        equal_min=inclusive_min)

    return ret


def get_info_value(variant: Dict[str, Any], key: str) -> Optional[float]:
    for item in variant['INFO'].split(';'):
        if '=' in item:
            k, v = item.split('=')
            if k == key:
                return float(v)
    return None
