from typing import List, IO
from .tools import edit_fpath
from .template import Processor
from .index_files import SamtoolsIndexFa, GATKCreateSequenceDictionary


class VariantFiltering(Processor):

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

        if self.variant_caller == 'mutect2':
            self.filter_mutect_calls()
            self.remove_variants()

        elif self.variant_caller == 'haplotype-caller':
            self.filter_haplotype_variants()
            self.remove_variants()

        else:
            self.logger.info(f"variant_caller = '{self.variant_caller}', skip variant filtering")

        return self.vcf

    def filter_mutect_calls(self):
        self.vcf = FilterMutectCalls(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.ref_fa)

    def filter_haplotype_variants(self):
        self.vcf = FilterHaplotypeVariants(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.ref_fa)

    def remove_variants(self):
        self.vcf = RemoveVariants(self.settings).main(
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
        log = f'{self.outdir}/gatk-FilterMutectCalls.log'

        self.call(self.CMD_LINEBREAK.join([
            'gatk FilterMutectCalls',
            f'--variant {self.vcf}',
            f'--reference {self.ref_fa}',
            f'--output {output}',
            f'--filtering-stats {stats_tsv}',
            f'--create-output-variant-index false',
            f'1> {log}',
            f'2> {log}',
        ]))

        self.vcf = output


class FilterHaplotypeVariants(Processor):

    vcf: str
    ref_fa: str

    snp_vcf: str
    indel_vcf: str
    flagged_snp_vcf: str
    flagged_indel_vcf: str
    flagged_combined_vcf: str

    def main(self, vcf: str, ref_fa: str) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa

        self.set_vcf_paths()
        self.separate_snp_indel()
        self.snp_variant_filtration()
        self.indel_variant_filtration()
        self.combine_flagged_snp_indel()

        return self.flagged_combined_vcf

    def set_vcf_paths(self):
        self.snp_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-snp.vcf',
            dstdir=self.workdir)

        self.indel_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-indel.vcf',
            dstdir=self.workdir)

        self.flagged_snp_vcf = edit_fpath(
            fpath=self.snp_vcf,
            old_suffix='.vcf',
            new_suffix='-flagged.vcf',
            dstdir=self.workdir)

        self.flagged_indel_vcf = edit_fpath(
            fpath=self.indel_vcf,
            old_suffix='.vcf',
            new_suffix='-flagged.vcf',
            dstdir=self.workdir)

        self.flagged_combined_vcf = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-snp-indel-flagged.vcf',
            dstdir=self.workdir)

    def separate_snp_indel(self):
        log = f'{self.outdir}/gatk-SelectVariants.log'
        self.call(self.CMD_LINEBREAK.join([
            'gatk SelectVariants',
            f'-variant {self.vcf}',
            f'--select-type-to-include SNP',
            f'-output {self.snp_vcf}',
            f'1>> {log}',
            f'2>> {log}',
        ]))
        self.call(self.CMD_LINEBREAK.join([
            'gatk SelectVariants',
            f'--variant {self.vcf}',
            f'--select-type-to-include INDEL',
            f'--output {self.indel_vcf}',
            f'1>> {log}',
            f'2>> {log}',
        ]))

    def snp_variant_filtration(self):
        log = f'{self.outdir}/gatk-VariantFiltration.log'
        self.call(self.CMD_LINEBREAK.join([
            'gatk VariantFiltration',
            f'--variant {self.snp_vcf}',
            f'-filter "QD < 2.0" --filter-name "QD2"',
            f'-filter "QD < 2.0" --filter-name "QD2"',
            f'-filter "QUAL < 30.0" --filter-name "QUAL30"',
            f'-filter "SOR > 3.0" --filter-name "SOR3"',
            f'-filter "FS > 60.0" --filter-name "FS60"',
            f'-filter "MQ < 40.0" --filter-name "MQ40"',
            f'-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5"',
            f'-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"',
            f'--output {self.flagged_snp_vcf}',
            f'1>> {log}',
            f'2>> {log}',
        ]))

    def indel_variant_filtration(self):
        log = f'{self.outdir}/gatk-VariantFiltration.log'
        self.call(self.CMD_LINEBREAK.join([
            'gatk VariantFiltration',
            f'--variant {self.indel_vcf}',
            f'-filter "QD < 2.0" --filter-name "QD2"',
            f'-filter "QUAL < 30.0" --filter-name "QUAL30"',
            f'-filter "FS > 200.0" --filter-name "FS200"',
            f'-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"',
            f'--output {self.flagged_indel_vcf}',
            f'1>> {log}',
            f'2>> {log}',
        ]))

    def combine_flagged_snp_indel(self):
        log = f'{self.outdir}/gatk-MergeVcfs.log'
        self.call(self.CMD_LINEBREAK.join([
            'gatk MergeVcfs',
            f'--INPUT {self.flagged_snp_vcf}',
            f'--INPUT {self.flagged_indel_vcf}',
            f'--OUTPUT {self.flagged_combined_vcf}',
            f'1>> {log}',
            f'2>> {log}',
        ]))


class RemoveVariants(Processor):

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
        self.write_to_output_vcf()
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

    def write_to_output_vcf(self):
        total, passed = 0, 0
        for line in self.reader:
            if line.startswith('#'):  # header
                self.writer.write(line)
                continue

            # variant line
            total += 1
            if self.__passed(line):
                passed += 1
                self.writer.write(line)

        self.__log_result(total, passed)

    def __passed(self, line: str) -> bool:
        variant_flags = line.split('\t')[6].split(';')
        red_flags = set(variant_flags).intersection(set(self.flags))
        return len(red_flags) == 0

    def __log_result(self, total: int, passed: int):
        percentage = passed / total * 100
        msg = f'''\
Remove variants having any one of the following flags: {', '.join(self.flags)}
Total variants: {total}
Remaining variants: {passed} ({percentage:.2f}%)'''
        self.logger.info(msg)

    def close_files(self):
        self.reader.close()
        self.writer.close()
