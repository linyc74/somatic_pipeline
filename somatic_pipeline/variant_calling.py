import os
from abc import ABC
from os.path import samefile, dirname
from typing import Optional, List, Dict, Callable, Tuple
from .tools import edit_fpath
from .template import Processor
from .vcf2csv import Vcf2Csv
from .constant import TUMOR, NORMAL
from .variant_filtering import FlagVariants, RemoveVariants
from .index_files import SamtoolsIndexFa, GATKIndexVcf, GATKCreateSequenceDictionary, SamtoolsIndexBam, BgzipIndex


class Params:

    def __init__(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            call_region_bed: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],
            variant_flagging_criteria: Optional[str],
            variant_removal_flags: List[str],
            only_pass: bool):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.call_region_bed = call_region_bed
        self.panel_of_normal_vcf = panel_of_normal_vcf
        self.germline_resource_vcf = germline_resource_vcf
        self.variant_flagging_criteria = variant_flagging_criteria
        self.variant_removal_flags = variant_removal_flags
        self.only_pass = only_pass


class VariantCalling(Processor):

    DSTDIR_NAME = 'callers'
    TN_PAIRED_MODE = 'tumor-normal paired'
    TUMOR_ONLY_MODE = 'tumor-only'

    variant_callers: List[str]
    params: Params

    dstdir: str
    mode: str
    mode_caller_to_method: Dict[str, Dict[str, Callable]]
    vcfs: List[str]

    def main(
            self,
            variant_callers: List[str],
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            call_region_bed: Optional[str],
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str],
            variant_flagging_criteria: Optional[str],
            variant_removal_flags: List[str],
            only_pass: bool) -> List[str]:

        self.variant_callers = variant_callers
        self.params = Params(
            ref_fa=ref_fa,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            call_region_bed=call_region_bed,
            panel_of_normal_vcf=panel_of_normal_vcf,
            germline_resource_vcf=germline_resource_vcf,
            variant_flagging_criteria=variant_flagging_criteria,
            variant_removal_flags=variant_removal_flags,
            only_pass=only_pass)

        self.index_ref_fa_and_bams()
        self.make_dstdir()

        self.vcfs = []
        for caller in self.variant_callers:

            method = self.get_caller_method(caller)
            flagged_vcf, filtered_vcf = method(self.params)

            self.move_and_process(vcf=flagged_vcf, dstname=caller)
            vcf_gz = self.move_and_process(vcf=filtered_vcf, dstname=f'{caller}-filtered')

            self.vcfs.append(vcf_gz)

        self.remove_bams()

        return self.vcfs

    def index_ref_fa_and_bams(self):
        SamtoolsIndexFa(self.settings).main(fa=self.params.ref_fa)
        SamtoolsIndexBam(self.settings).main(bam=self.params.tumor_bam)
        if self.params.normal_bam is not None:
            SamtoolsIndexBam(self.settings).main(bam=self.params.normal_bam)

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def get_caller_method(self, caller: str) -> Callable[[Params], Tuple[str, str]]:
        if self.params.normal_bam is not None:
            mode = 'tumor-normal paired'
            ret = {
                'mutect2': Mutect2TNPaired(self.settings).main,
                'muse': Muse(self.settings).main,
                'varscan': Varscan(self.settings).main,
                'vardict': VarDictTNPaired(self.settings).main,
                'lofreq': LoFreqTNPaired(self.settings).main,
                'somatic-sniper': SomaticSniper(self.settings).main,
            }.get(caller, None)
        else:
            mode = 'tumor-only'
            ret = {
                'mutect2': Mutect2TumorOnly(self.settings).main,
                'haplotype-caller': HaplotypeCaller(self.settings).main,
                'vardict': VarDictTumorOnly(self.settings).main,
                'lofreq': LoFreqTumorOnly(self.settings).main,
            }.get(caller, None)

        assert ret is not None, f'For {mode} mode, variant caller "{caller}" does not work'

        return ret

    def move_and_process(self, vcf: str, dstname: str) -> str:
        dst_vcf = f'{self.dstdir}/{dstname}.vcf'
        self.call(f'mv {vcf} {dst_vcf}')
        Vcf2Csv(self.settings).main(vcf=dst_vcf)
        vcf_gz = BgzipIndex(self.settings).main(vcf=dst_vcf, keep=False)
        return vcf_gz

    def remove_bams(self):
        for bam in [self.params.tumor_bam, self.params.normal_bam]:
            if bam is None:
                continue
            if samefile(dirname(bam), self.workdir):  # bam in the workdir
                self.call(f'rm {bam}')


#


class CallerBase(Processor, ABC):

    params: Params

    vcf: str  # dynamically track the VCF file path
    flagged_vcf: str  # static, the VCF file before variant removal
    filtered_vcf: str  # static, the VCF file after variant removal

    def main(self, params: Params) -> Tuple[str, str]:
        self.params = params
        self.main__()
        return self.flagged_vcf, self.filtered_vcf

    def main__(self):
        pass

    def flag_and_remove_variants(self):
        if self.params.variant_flagging_criteria is not None:
            self.vcf = FlagVariants(self.settings).main(
                vcf=self.vcf,
                variant_flagging_criteria=self.params.variant_flagging_criteria)

        self.flagged_vcf = self.vcf

        self.vcf = RemoveVariants(self.settings).main(
            vcf=self.vcf,
            flags=self.params.variant_removal_flags,
            only_pass=self.params.only_pass)

        self.filtered_vcf = self.vcf


#


class GATKBase(CallerBase):

    def create_sequence_dictionary(self):
        GATKCreateSequenceDictionary(self.settings).main(ref_fa=self.params.ref_fa)


class Mutect2Base(GATKBase):

    pon_args: List[str]
    germline_resource_args: List[str]
    f1r2_tar_gz: str

    def prepare_mutect2_resource_vcfs(self):
        self.pon_args, self.germline_resource_args = PrepareMutect2ResourceVcfs(self.settings).main(
            panel_of_normal_vcf=self.params.panel_of_normal_vcf,
            germline_resource_vcf=self.params.germline_resource_vcf)

    def filter_mutect_calls(self):
        self.vcf = FilterMutectCalls(self.settings).main(
            vcf=self.vcf,
            ref_fa=self.params.ref_fa,
            f1r2_tar_gz=self.f1r2_tar_gz)


class PrepareMutect2ResourceVcfs(Processor):

    panel_of_normal_vcf: Optional[str]
    germline_resource_vcf: Optional[str]

    pon_args: List[str]
    germline_resource_args: List[str]

    def main(
            self,
            panel_of_normal_vcf: Optional[str],
            germline_resource_vcf: Optional[str]) -> Tuple[List[str], List[str]]:

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

    def main__(self):
        self.create_sequence_dictionary()
        self.prepare_mutect2_resource_vcfs()
        self.mutect2()
        self.filter_mutect_calls()
        self.flag_and_remove_variants()

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        bed = [] if self.params.call_region_bed is None else [f'--intervals {self.params.call_region_bed}']
        self.vcf = f'{self.workdir}/mutect2.vcf'
        self.f1r2_tar_gz = f'{self.workdir}/gatk-mutect2-f1r2.tar.gz'

        args = [
            'gatk Mutect2',
            f'--reference {self.params.ref_fa}',
            f'--input {self.params.tumor_bam}',
            f'--input {self.params.normal_bam}',
            f'--tumor-sample {TUMOR}',
            f'--normal-sample {NORMAL}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'--f1r2-tar-gz {self.f1r2_tar_gz}',
        ] + bed + self.pon_args + self.germline_resource_args + [
            f'1> {log}',
            f'2> {log}',
        ]

        self.call(self.CMD_LINEBREAK.join(args))


class Mutect2TumorOnly(Mutect2Base):

    def main__(self):
        self.create_sequence_dictionary()
        self.prepare_mutect2_resource_vcfs()
        self.mutect2()
        self.filter_mutect_calls()
        self.flag_and_remove_variants()

    def mutect2(self):
        log = f'{self.outdir}/gatk-Mutect2.log'
        bed = [] if self.params.call_region_bed is None else [f'--intervals {self.params.call_region_bed}']
        self.vcf = f'{self.workdir}/mutect2.vcf'
        self.f1r2_tar_gz = f'{self.workdir}/gatk-mutect2-f1r2.tar.gz'
        cmd = self.CMD_LINEBREAK.join([
            'gatk Mutect2',
            f'--reference {self.params.ref_fa}',
            f'--input {self.params.tumor_bam}',
            f'--output {self.vcf}',
            f'--native-pair-hmm-threads {self.threads}',
            f'--f1r2-tar-gz {self.f1r2_tar_gz}',
        ] + bed + self.pon_args + self.germline_resource_args + [
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class HaplotypeCaller(GATKBase):

    def main__(self):
        self.create_sequence_dictionary()
        self.haplotype_caller()
        self.filter_haplotype_variants()
        self.flag_and_remove_variants()

    def haplotype_caller(self):
        log = f'{self.outdir}/gatk-HaplotypeCaller.log'
        bed = [] if self.params.call_region_bed is None else [f'--intervals {self.params.call_region_bed}']
        output = f'{self.workdir}/haplotype-caller.vcf'
        args = [
            'gatk HaplotypeCaller',
            f'--reference {self.params.ref_fa}',
            f'--input {self.params.tumor_bam}',
        ] + bed + [
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
            ref_fa=self.params.ref_fa)


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


class Muse(CallerBase):

    WES_OR_WGS = '-E'  # '-E' for WES; '-G' for WGS

    call_result_txt: str

    def main__(self):
        self.clean_up_call_region_bed()
        self.muse_call()
        self.muse_sump()
        self.flag_and_remove_variants()

    def clean_up_call_region_bed(self):
        if self.params.call_region_bed is not None:
            self.params.call_region_bed = CleanUpCallRegionBed(self.settings).main(
                bed=self.params.call_region_bed)

    def muse_call(self):
        output = f'{self.workdir}/MuSE-call-result'
        bed = [] if self.params.call_region_bed is None else [f'-l {self.params.call_region_bed}']
        log = f'{self.outdir}/MuSE-call.log'
        cmd = self.CMD_LINEBREAK.join([
            'MuSE call',
            f'-f {self.params.ref_fa}',
            f'-O {output}',
        ] + bed + [
            self.params.tumor_bam,
            self.params.normal_bam,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)
        self.call_result_txt = output + '.MuSE.txt'

    def muse_sump(self):
        self.vcf = f'{self.workdir}/muse.vcf'
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


class CleanUpCallRegionBed(Processor):

    bed: str
    output_bed: str

    def main(self, bed: str):
        self.bed = bed

        self.output_bed = f'{self.workdir}/call-region.bed'

        self.logger.info(f'Cleaning up call region BED file by only using the first 3 fields:\n"{self.bed}" --> "{self.output_bed}"')

        with open(self.bed) as reader:
            with open(self.output_bed, 'w') as writer:
                for line in reader:
                    items = line.strip().split('\t')
                    new_line = '\t'.join(items[0:3])
                    writer.write(new_line + '\n')

        return self.output_bed


#


class Varscan(CallerBase):

    normal_pileup: str
    tumor_pileup: str
    snp_vcf: str
    indel_vcf: str

    def main__(self):
        self.samtools_mpileup()
        self.varscan_somatic()
        self.compress_vcfs()
        self.index_vcfs()
        self.concat_snp_indel_vcfs()
        self.flag_and_remove_variants()

    def samtools_mpileup(self):
        self.tumor_pileup = f'{self.workdir}/{TUMOR}-pileup'
        self.normal_pileup = f'{self.workdir}/{NORMAL}-pileup'
        for bam, output, name in [
            (self.params.tumor_bam, self.tumor_pileup, TUMOR),
            (self.params.normal_bam, self.normal_pileup, NORMAL),
        ]:
            self.__mpileup(bam=bam, output=output, name=name)

    def __mpileup(self, bam: str, output: str, name: str):
        log = f'{self.outdir}/samtools-mpileup-{name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'samtools mpileup',
            f'--fasta-ref {self.params.ref_fa}',
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
        self.vcf = f'{self.workdir}/varscan.vcf'
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


class VarDictBase(CallerBase):

    ALLELE_FREQ_THRESHOLD = 0.01

    genome_file: str

    def build_call_region_bed_for_wgs_mode(self):
        self.__write_genome_file()
        self.__write_bed()

    def __write_genome_file(self):
        chr_to_length = {}
        this_chr = None
        with open(self.params.ref_fa) as reader:
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

    def __write_bed(self):
        self.params.call_region_bed = f'{self.workdir}/vardict-wgs.bed'
        # Create regions with a window size of 50150 bp and overlapping size 150 bp
        self.call(f'bedtools makewindows -g {self.genome_file} -w 50150 -s 50000 > {self.params.call_region_bed}')


class VarDictTumorOnly(VarDictBase):

    def main__(self):
        if self.params.call_region_bed is None:
            self.build_call_region_bed_for_wgs_mode()
        self.run_vardict()
        self.flag_and_remove_variants()

    def run_vardict(self):
        log = f'{self.outdir}/vardict.log'
        self.vcf = f'{self.workdir}/vardict.vcf'
        args = [
            'VarDict',
            f'-G {self.params.ref_fa}',
            f'-f {self.ALLELE_FREQ_THRESHOLD}',
            f'-N {TUMOR}',
            f'-b {self.params.tumor_bam}',
            f'-th {self.threads}',
            '-c 1',  # column for chromosome
            '-S 2',  # column for region start
            '-E 3',  # column for region end
            '-g 4',  # column for gene name
            self.params.call_region_bed,
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

    def main__(self):
        if self.params.call_region_bed is None:
            self.build_call_region_bed_for_wgs_mode()
        self.run_vardict()
        self.flag_and_remove_variants()

    def run_vardict(self):
        log = f'{self.outdir}/vardict.log'
        self.vcf = f'{self.workdir}/vardict.vcf'
        args = [
            'VarDict',
            f'-G {self.params.ref_fa}',
            f'-f {self.ALLELE_FREQ_THRESHOLD}',
            f'-N {TUMOR}',
            f'-b "{self.params.tumor_bam}|{self.params.normal_bam}"',
            f'-th {self.threads}',
            '-c 1',  # column for chromosome
            '-S 2',  # column for region start
            '-E 3',  # column for region end
            '-g 4',  # column for gene name
            self.params.call_region_bed,
            f'2> {log}',
            '|',
            'testsomatic.R',
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


class LoFreqTumorOnly(CallerBase):

    def main__(self):
        self.lofreq_call_parallel()
        self.flag_and_remove_variants()

    def lofreq_call_parallel(self):
        self.vcf = f'{self.workdir}/lofreq.vcf'
        log = f'{self.outdir}/lofreq-call-parallel.log'
        bed = [] if self.params.call_region_bed is None else [f'--bed {self.params.call_region_bed}']
        args = [
            'lofreq call-parallel',
            f'--ref {self.params.ref_fa}',
        ] + bed + [
            f'--pp-threads {self.threads}',
            '--call-indels',  # works only when BAM file contains indel qualities, GATK BQSR does that
            f'--out {self.vcf}',
            self.params.tumor_bam,
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))


class LoFreqTNPaired(CallerBase):

    def main__(self):
        self.lofreq_somatic()
        self.concat_snp_indel_vcfs()
        self.flag_and_remove_variants()

    def lofreq_somatic(self):
        log = f'{self.outdir}/lofreq-somatic.log'
        bed = [] if self.params.call_region_bed is None else [f'--bed {self.params.call_region_bed}']
        args = [
            'lofreq somatic',
            f'--normal {self.params.normal_bam}',
            f'--tumor {self.params.tumor_bam}',
            f'--ref {self.params.ref_fa}',
        ] + bed + [
            f'--threads {self.threads}',
            '--call-indels',  # works only when BAM file contains indel qualities, GATK BQSR does that
            f'-o {self.workdir}/lofreq-',
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))

    def concat_snp_indel_vcfs(self):
        snp_vcf = f'{self.workdir}/lofreq-somatic_final.snvs.vcf.gz'
        indel_vcf = f'{self.workdir}/lofreq-somatic_final.indels.vcf.gz'
        self.vcf = f'{self.workdir}/lofreq.vcf'
        log = f'{self.outdir}/bcftools-merge.log'
        cmd = self.CMD_LINEBREAK.join([
            'bcftools concat',
            '--allow-overlaps',  # first coordinate of the next file can precede last record of the current file
            f'--threads {self.threads}',
            f'--output-type v',  # uncompressed VCF
            f'--output {self.vcf}',
            snp_vcf,
            indel_vcf,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


#


class SomaticSniper(CallerBase):

    def main__(self):
        self.run_somatic_sniper()
        self.flag_and_remove_variants()

    def run_somatic_sniper(self):
        self.vcf = f'{self.workdir}/somatic-sniper.vcf'
        log = f'{self.outdir}/bam-somaticsniper.log'
        args = [
            'bam-somaticsniper',
            f'-n {NORMAL}',
            f'-t {TUMOR}',
            f'-F vcf',  # output format
            f'-f {self.params.ref_fa}',
            self.params.tumor_bam,
            self.params.normal_bam,
            self.vcf,
            f'1> {log}',
            f'2> {log}',
        ]
        self.call(self.CMD_LINEBREAK.join(args))
