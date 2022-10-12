from abc import ABC
from typing import Optional, List, IO
from .tools import edit_fpath
from .template import Processor
from .constant import TUMOR, NORMAL
from .index_files import SamtoolsIndexFa, GATKIndexVcf, GATKCreateSequenceDictionary, SamtoolsIndexBam


class VariantCalling(Processor):

    variant_caller: str
    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]
    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]
    vardict_call_region_bed: Optional[str]
    variant_removal_flags: List[str]

    vcf: str

    def main(
            self,
            variant_caller: str,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],
            vardict_call_region_bed: Optional[str],
            variant_removal_flags: List[str]) -> str:

        self.variant_caller = variant_caller
        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf
        self.vardict_call_region_bed = vardict_call_region_bed
        self.variant_removal_flags = variant_removal_flags

        tumor_normal_paired = normal_bam is not None

        if tumor_normal_paired:
            if variant_caller == 'mutect2':
                self.mutect2_tn_paired()
            elif variant_caller == 'muse':
                self.muse()
            elif variant_caller == 'varscan':
                self.varscan()
            elif variant_caller == 'vardict':
                self.vardict_tn_paired()
            else:
                raise AssertionError(f'Variant caller "{variant_caller}" not available for tumor-normal paired mode')

        else:  # tumor-only
            if variant_caller == 'mutect2':
                self.mutect2_tumor_only()
            elif variant_caller == 'haplotype-caller':
                self.haplotype_caller()
            elif variant_caller == 'vardict':
                self.vardict_tumor_only()
            else:
                raise AssertionError(f'Variant caller "{variant_caller}" not available for tumor-only mode')

        return self.vcf

    def mutect2_tn_paired(self):
        self.vcf = Mutect2TNPaired(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            germline_resource_vcf=self.germline_resource_vcf,
            variant_removal_flags=self.variant_removal_flags)

    def mutect2_tumor_only(self):
        self.vcf = Mutect2TumorOnly(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            germline_resource_vcf=self.germline_resource_vcf,
            variant_removal_flags=self.variant_removal_flags)

    def haplotype_caller(self):
        self.vcf = HaplotypeCaller(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            variant_removal_flags=self.variant_removal_flags)

    def muse(self):
        self.vcf = Muse(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            variant_removal_flags=self.variant_removal_flags)

    def varscan(self):
        self.vcf = Varscan(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            variant_removal_flags=self.variant_removal_flags)

    def vardict_tn_paired(self):
        self.vcf = VarDictTNPaired(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            bed=self.vardict_call_region_bed,
            variant_removal_flags=self.variant_removal_flags)

    def vardict_tumor_only(self):
        self.vcf = VarDictTumorOnly(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            bed=self.vardict_call_region_bed,
            variant_removal_flags=self.variant_removal_flags)


class Base(Processor, ABC):

    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]
    variant_removal_flags: List[str]

    vcf: str

    def index_ref_fa_and_bams(self):
        SamtoolsIndexFa(self.settings).main(fa=self.ref_fa)
        SamtoolsIndexBam(self.settings).main(bam=self.tumor_bam)
        if self.normal_bam is not None:
            SamtoolsIndexBam(self.settings).main(bam=self.normal_bam)

    def remove_variants(self):
        if len(self.variant_removal_flags) > 0:
            self.vcf = RemoveVariants(self.settings).main(
                vcf=self.vcf,
                flags=self.variant_removal_flags)


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


#


class GATKBase(Base):

    def create_sequence_dictionary(self):
        GATKCreateSequenceDictionary(self.settings).main(ref_fa=self.ref_fa)


class Mutect2Base(GATKBase):

    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    pon_args: List[str]
    germline_resource_args: List[str]
    f1r2_tar_gz: str

    def prepare_mutect2_resource_vcfs(self):
        self.pon_args, self.germline_resource_args = PrepareMutect2ResourceVcfs(self.settings).main(
            panel_of_normal_vcf=self.panel_of_normal_vcf,
            germline_resource_vcf=self.germline_resource_vcf)

    def filter_mutect_calls(self):
        self.vcf = FilterMutectCalls(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.ref_fa,
            f1r2_tar_gz=self.f1r2_tar_gz)


class PrepareMutect2ResourceVcfs(Processor):

    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    pon_args: List[str]
    germline_resource_args: List[str]

    def main(
            self,
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str]):

        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf

        self.prepare_panel_of_normal()
        self.prepare_germline_resource()

        return self.pon_args, self.germline_resource_args

    def prepare_panel_of_normal(self):
        if self.panel_of_normal_vcf is None:
            self.pon_args = []
        else:
            self.panel_of_normal_vcf = self.__copy(self.panel_of_normal_vcf)
            self.__index(self.panel_of_normal_vcf)
            self.pon_args = [f'--panel-of-normals {self.panel_of_normal_vcf}']

    def prepare_germline_resource(self):
        if self.germline_resource_vcf is None:
            self.germline_resource_args = []
        else:
            self.germline_resource_vcf = self.__copy(self.germline_resource_vcf)
            self.__index(self.germline_resource_vcf)
            self.germline_resource_args = [f'--germline-resource {self.germline_resource_vcf}']

    def __copy(self, vcf: str) -> str:
        dst = edit_fpath(fpath=vcf, dstdir=self.workdir)
        self.call(f'cp {vcf} {dst}')
        return dst

    def __index(self, vcf: str):
        GATKIndexVcf(self.settings).main(vcf=vcf)


class FilterMutectCalls(Processor):

    vcf: str
    ref_fa: str
    f1r2_tar_gz: str

    orientation_artifact_tar_gz: str

    def main(self, vcf: str, ref_fa: str, f1r2_tar_gz: str) -> str:

        self.vcf = vcf
        self.ref_fa = ref_fa
        self.f1r2_tar_gz = f1r2_tar_gz

        self.learn_read_orientation_model()
        self.filter_mutect_calls()

        return self.vcf

    def learn_read_orientation_model(self):
        log = f'{self.outdir}/gatk-LearnReadOrientationModel.log'
        self.orientation_artifact_tar_gz = f'{self.workdir}/artifact-prior-table.tar.gz'
        self.call(self.CMD_LINEBREAK.join([
            'gatk LearnReadOrientationModel',
            f'--input {self.f1r2_tar_gz}',
            f'--output {self.orientation_artifact_tar_gz}',
            f'1> {log}',
            f'2> {log}',
        ]))

    def filter_mutect_calls(self):
        log = f'{self.outdir}/gatk-FilterMutectCalls.log'
        stats_tsv = f'{self.outdir}/filter-mutect-stats.tsv'
        output = edit_fpath(
            fpath=self.vcf,
            old_suffix='.vcf',
            new_suffix='-filter-mutect-calls.vcf',
            dstdir=self.workdir)
        self.call(self.CMD_LINEBREAK.join([
            'gatk FilterMutectCalls',
            f'--variant {self.vcf}',
            f'--reference {self.ref_fa}',
            f'--output {output}',
            f'--filtering-stats {stats_tsv}',
            f'--orientation-bias-artifact-priors {self.orientation_artifact_tar_gz}',
            f'--create-output-variant-index false',
            f'1> {log}',
            f'2> {log}',
        ]))
        self.vcf = output


class Mutect2TNPaired(Mutect2Base):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()
        self.create_sequence_dictionary()
        self.prepare_mutect2_resource_vcfs()
        self.mutect2()
        self.filter_mutect_calls()
        self.remove_variants()

        return self.vcf

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        self.vcf = f'{self.workdir}/raw.vcf'
        self.f1r2_tar_gz = f'{self.workdir}/gatk-mutect2-f1r2.tar.gz'
        args = [
            'gatk Mutect2',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_bam}',
            f'--input {self.normal_bam}',
            f'--tumor-sample {TUMOR}',
            f'--normal-sample {NORMAL}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'--f1r2-tar-gz {self.f1r2_tar_gz}',
        ] + self.pon_args + self.germline_resource_args + [
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))


class Mutect2TumorOnly(Mutect2Base):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = None
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()
        self.create_sequence_dictionary()
        self.prepare_mutect2_resource_vcfs()
        self.mutect2()
        self.filter_mutect_calls()
        self.remove_variants()

        return self.vcf

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        self.vcf = f'{self.workdir}/raw.vcf'
        self.f1r2_tar_gz = f'{self.workdir}/gatk-mutect2-f1r2.tar.gz'
        cmd = self.CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_bam}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'--f1r2-tar-gz {self.f1r2_tar_gz}',
        ] + self.pon_args + self.germline_resource_args + [
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class HaplotypeCaller(GATKBase):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()
        self.create_sequence_dictionary()
        self.haplotype_caller()
        self.filter_haplotype_variants()
        self.remove_variants()

        return self.vcf

    def haplotype_caller(self):
        log = f'{self.outdir}/gatk-HaplotypeCaller.log'
        output = f'{self.workdir}/raw.vcf'
        args = [
            'gatk HaplotypeCaller',
            f'--reference {self.ref_fa}',
            f'--input {self.tumor_bam}',
            f'--output {output}',
            f'--native-pair-hmm-threads {self.threads}',
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))
        self.vcf = output

    def filter_haplotype_variants(self):
        self.vcf = FilterHaplotypeVariants(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.ref_fa)


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


#


class Muse(Base):

    WES_OR_WGS = '-E'  # '-E' for WES; '-G' for WGS

    call_result_txt: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()
        self.muse_call()
        self.muse_sump()
        self.remove_variants()

        return self.vcf

    def muse_call(self):
        log = f'{self.outdir}/MuSE-call.log'
        output = f'{self.workdir}/MuSE-call-result'
        cmd = self.CMD_LINEBREAK.join([
            'MuSE call',
            f'-f {self.ref_fa}',
            f'-O {output}',
            self.tumor_bam,
            self.normal_bam,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)
        self.call_result_txt = output + '.MuSE.txt'

    def muse_sump(self):
        self.vcf = f'{self.workdir}/raw.vcf'
        log = f'{self.outdir}/MuSE-sump.log'
        cmd = self.CMD_LINEBREAK.join([
            'MuSE sump',
            f'-I {self.call_result_txt}',
            self.WES_OR_WGS,
            f'-O {self.vcf}',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Varscan(Base):

    normal_pileup: str
    tumor_pileup: str
    snp_vcf: str
    indel_vcf: str

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()
        self.samtools_mpileup()
        self.varscan_somatic()
        self.compress_vcfs()
        self.index_vcfs()
        self.concat_snp_indel_vcfs()
        self.remove_variants()

        return self.vcf

    def samtools_mpileup(self):
        self.tumor_pileup = f'{self.workdir}/{TUMOR}-pileup'
        self.normal_pileup = f'{self.workdir}/{NORMAL}-pileup'
        for bam, output, name in [
            (self.tumor_bam, self.tumor_pileup, TUMOR),
            (self.normal_bam, self.normal_pileup, NORMAL),
        ]:
            self.__mpileup(bam=bam, output=output, name=name)

    def __mpileup(self, bam: str, output: str, name: str):
        log = f'{self.outdir}/samtools-mpileup-{name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'samtools mpileup',
            f'--fasta-ref {self.ref_fa}',
            f'--output {output}',
            bam,
            f'2> {log}',
        ])
        self.call(cmd)

    def varscan_somatic(self):
        self.snp_vcf = f'{self.workdir}/varscan-snp.vcf'
        self.indel_vcf = f'{self.workdir}/varscan-indel.vcf'
        log = f'{self.outdir}/varscan-somatic.log'
        cmd = self.CMD_LINEBREAK.join([
            'varscan somatic',
            self.normal_pileup,
            self.tumor_pileup,
            f'--output-snp {self.snp_vcf}',
            f'--output-indel {self.indel_vcf}',
            '--strand-filter 1',  # removes variants with >90% strand bias
            '--output-vcf 1',  # output VCF instead of VarScan native format
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)

    def compress_vcfs(self):
        for vcf in [self.snp_vcf, self.indel_vcf]:
            self.call(f'bgzip --force --threads {self.threads} {vcf}')
        self.snp_vcf += '.gz'
        self.indel_vcf += '.gz'

    def index_vcfs(self):
        for vcf in [self.snp_vcf, self.indel_vcf]:
            self.call(f'bcftools index --force --threads {self.threads} {vcf}')

    def concat_snp_indel_vcfs(self):
        log = f'{self.outdir}/bcftools-merge.log'
        self.vcf = f'{self.workdir}/raw.vcf'
        cmd = self.CMD_LINEBREAK.join([
            'bcftools concat',
            '--allow-overlaps',  # first coordinate of the next file can precede last record of the current file
            f'--threads {self.threads}',
            f'--output-type v',  # uncompressed VCF
            f'--output {self.vcf}',
            self.snp_vcf,
            self.indel_vcf,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


#


class VarDictBase(Base):

    ALLELE_FREQ_THRESHOLD = 0.01

    bed: str
    genome_file: str

    def build_bed_for_wgs_mode(self):
        self.__genome_file()
        self.__bed()

    def __genome_file(self):
        chr_to_length = {}
        this_chr = None
        with open(self.ref_fa) as reader:
            for line in reader:
                if line.startswith('>'):
                    this_chr = line[1:].split(' ')[0]
                    chr_to_length.setdefault(this_chr, 0)
                else:
                    chr_to_length[this_chr] += len(line.strip())

        self.genome_file = f'{self.workdir}/vardict-genome-file.txt'
        with open(self.genome_file, 'w') as writer:
            for c, l in chr_to_length.items():
                writer.write(f'{c}\t{l}\n')

    def __bed(self):
        self.bed = f'{self.workdir}/vardict-wgs.bed'
        # Create regions with a window size of 50150 bp and overlapping size 150 bp
        self.call(f'bedtools makewindows -g {self.genome_file} -w 50150 -s 50000 > {self.bed}')


class VarDictTumorOnly(VarDictBase):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            variant_removal_flags: List[str],
            bed: Optional[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = None
        self.variant_removal_flags = variant_removal_flags
        self.bed = bed

        self.index_ref_fa_and_bams()
        if self.bed is None:
            self.build_bed_for_wgs_mode()
        self.run_vardict()
        self.remove_variants()

        return self.vcf

    def run_vardict(self):
        log = f'{self.outdir}/vardict.log'
        self.vcf = f'{self.workdir}/raw.vcf'
        args = [
            'vardict',
            f'-G {self.ref_fa}',
            f'-f {self.ALLELE_FREQ_THRESHOLD}',
            f'-N {TUMOR}',
            f'-b {self.tumor_bam}',
            '-c 1',  # column for chromosome
            '-S 2',  # column for region start
            '-E 3',  # column for region end
            '-g 4',  # column for gene name
            self.bed,
            f'2> {log}',
            '|',
            'teststrandbias.R',
            f'2> {log}',
            '|',
            'var2vcf_valid.pl',
            f'-N {TUMOR}',
            '-E',  # do not print END tag
            f'-f {self.ALLELE_FREQ_THRESHOLD}',
            f'1> {self.vcf}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))


class VarDictTNPaired(VarDictBase):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            variant_removal_flags: List[str],
            bed: Optional[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags
        self.bed = bed

        self.index_ref_fa_and_bams()
        if self.bed is None:
            self.build_bed_for_wgs_mode()
        self.run_vardict()
        self.remove_variants()

        return self.vcf

    def run_vardict(self):
        log = f'{self.outdir}/vardict.log'
        self.vcf = f'{self.workdir}/raw.vcf'
        args = [
            'vardict',
            f'-G {self.ref_fa}',
            f'-f {self.ALLELE_FREQ_THRESHOLD}',
            f'-N {TUMOR}',
            f'-b "{self.tumor_bam}|{self.normal_bam}"',
            '-c 1',  # column for chromosome
            '-S 2',  # column for region start
            '-E 3',  # column for region end
            '-g 4',  # column for gene name
            self.bed,
            f'2> {log}',
            '|',
            'teststrandbias.R',
            f'2> {log}',
            '|',
            'var2vcf_paired.pl',
            f'-N "{TUMOR}|{NORMAL}"',
            f'-f {self.ALLELE_FREQ_THRESHOLD}',
            f'1> {self.vcf}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))


#


class LoFreq(Base):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()

        self.remove_variants()

        return self.vcf


class Platypus(Base):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()

        self.remove_variants()

        return self.vcf


class SomaticSniper(Base):

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            variant_removal_flags: List[str]) -> str:

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.variant_removal_flags = variant_removal_flags

        self.index_ref_fa_and_bams()

        self.remove_variants()

        return self.vcf
