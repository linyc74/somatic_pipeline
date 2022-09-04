from typing import List, IO
from .tools import edit_fpath
from .template import Processor
from .index_files import SamtoolsIndexFa, GATKCreateSequenceDictionary


class Mutect2VariantFiltering(Processor):

    vcf: str
    ref_fa: str
    variant_caller: str
    variant_removal_flags: List[str]

    def main(
            self,
            vcf: str,
            ref_fa: str,
            variant_caller: str,
            variant_removal_flags: List[str]) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa
        self.variant_caller = variant_caller
        self.variant_removal_flags = variant_removal_flags

        if self.valid_caller():
            self.filter_mutect_calls()
            self.variant_removal()
        else:
            self.logger.info(f"variant_caller = '{self.variant_caller}', skip mutect2 variant filtering")

        return self.vcf

    def valid_caller(self) -> bool:
        return self.variant_caller == 'mutect2'

    def filter_mutect_calls(self):
        self.vcf = FilterMutectCalls(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.ref_fa)

    def variant_removal(self):
        self.vcf = VariantRemoval(self.settings).main(
            vcf=self.vcf,
            flags=self.variant_removal_flags)


class FilterMutectCalls(Processor):

    vcf: str
    ref_fa: str

    def main(self, vcf: str, ref_fa: str) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa

        self.prepare_ref_fa()
        self.execute()

        return self.vcf

    def prepare_ref_fa(self):
        SamtoolsIndexFa(self.settings).main(fa=self.ref_fa)
        GATKCreateSequenceDictionary(self.settings).main(ref_fa=self.ref_fa)

    def execute(self):
        output = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-filter-mutect-calls.vcf',
            dstdir=self.workdir)

        stats_tsv = f'{self.outdir}/filter-mutect-stats.tsv'

        args = [
            'gatk FilterMutectCalls',
            f'--variant {self.vcf}',
            f'--reference {self.ref_fa}',
            f'--output {output}',
            f'--filtering-stats {stats_tsv}',
            f'--create-output-variant-index false',
        ]
        self.call(self.CMD_LINEBREAK.join(args))

        self.vcf = output


class VariantRemoval(Processor):

    vcf: str
    flags: List[str]

    output_vcf: str
    reader: IO
    writer: IO

    def main(self, vcf: str, flags: List[str]) -> str:

        self.vcf = vcf
        self.flags = flags

        self.set_output_vcf()
        self.open_files()

        for line in self.reader:
            if line.startswith('#'):  # header
                self.writer.write(line)
                continue
            if self.passed(line):
                self.writer.write(line)

        self.close_files()

        return self.output_vcf

    def set_output_vcf(self):
        self.output_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-variant-removal.vcf',
            dstdir=self.workdir)

    def open_files(self):
        self.reader = open(self.vcf)
        self.writer = open(self.output_vcf, 'w')

    def passed(self, line: str) -> bool:
        variant_flags = line.split('\t')[6].split(';')
        red_flags = set(variant_flags).intersection(set(self.flags))
        return len(red_flags) == 0

    def close_files(self):
        self.reader.close()
        self.writer.close()
