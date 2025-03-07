import argparse
from typing import List
from somatic_pipeline import Run


__VERSION__ = '1.12.0-beta'


PURPLE = '\033[95m'
CYAN = '\033[96m'
DARKCYAN = '\033[36m'
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
END = '\033[0m'


PROG = 'python somatic_pipeline'
DESCRIPTION = f'{BOLD}Somatic pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw){END}'
HELP_ARG = {
    'keys': ['-h', '--help'],
    'properties': {
        'action': 'help',
        'help': 'show this help message',
    }
}
VERSION_ARG = {
    'keys': ['-v', '--version'],
    'properties': {
        'action': 'version',
        'version': __VERSION__,
        'help': 'show version',
    }
}
DEBUG_ARG = {
    'keys': ['-d', '--debug'],
    'properties': {
        'action': 'store_true',
        'help': 'debug mode',
    }
}
GROUP_NAME_TO_ARGUMENTS = {
    'main':
        [
            {
                'keys': ['-r', '--ref-fa'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the reference genome fasta(.gz) file',
                }
            },
            {
                'keys': ['-1', '--tumor-fq1'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the tumor read 1 fastq(.gz) file',
                }
            },
            {
                'keys': ['-2', '--tumor-fq2'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the tumor read 2 fastq(.gz) file',
                }
            },
        ],

    'annotate':
        [
            {
                'keys': ['-f', '--vcf'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the vcf(.gz) file',
                }
            },
            {
                'keys': ['-r', '--ref-fa'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the reference genome fasta(.gz) file',
                }
            },
        ],

    'normal':
        [
            {
                'keys': ['-3', '--normal-fq1'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'path to the normal read 1 fastq(.gz) file (default: %(default)s)',
                }
            },
            {
                'keys': ['-4', '--normal-fq2'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'path to the normal read 2 fastq(.gz) file (default: %(default)s)',
                }
            },
        ],

    'general':
        [
            {
                'keys': ['-o', '--outdir'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'somatic_pipeline_outdir',
                    'help': 'path to the output directory (default: %(default)s)',
                }
            },
            {
                'keys': ['-t', '--threads'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 4,
                    'help': 'number of CPU threads (default: %(default)s)',
                }
            },
            DEBUG_ARG,
            HELP_ARG,
            VERSION_ARG,
        ],

    'pre-processing':
        [
            {
                'keys': ['--umi-length'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 0,
                    'help': 'UMI length (bp) to be removed (default: %(default)s)',
                }
            },
            {
                'keys': ['--clip-r1-5-prime'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 0,
                    'help': 'hard clip <int> bp from 5\' end of read 1 (default: %(default)s)',
                }
            },
            {
                'keys': ['--clip-r2-5-prime'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 0,
                    'help': 'hard clip <int> bp from 5\' end of read 2 (default: %(default)s)',
                }
            },
            {
                'keys': ['--read-aligner'],
                'properties': {
                    'type': str,
                    'required': False,
                    'choices': ['bwa', 'bowtie2'],
                    'default': 'bwa',
                    'help': 'read aligner (default: %(default)s)',
                }
            },
            {
                'keys': ['--skip-mark-duplicates'],
                'properties': {
                    'action': 'store_true',
                    'help': 'do not mark PCR duplicates',
                }
            },
            {
                'keys': ['--bqsr-known-variant-vcf'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'known variants VCF file for BQSR, if "None" then skip BQSR (default: %(default)s)',
                }
            },
            {
                'keys': ['--discard-bam'],
                'properties': {
                    'action': 'store_true',
                    'help': 'do not save sorted BAM files in outdir',
                }
            },
        ],

    'variant calling':
        [
            {
                'keys': ['--variant-callers'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'mutect2',
                    'help': 'comma-separated variant callers: mutect2,haplotype-caller,vardict,lofreq,muse,varscan,somatic-sniper, the former four works for tumor-only (default: %(default)s)',

                }
            },
            {
                'keys': ['--skip-variant-calling'],
                'properties': {
                    'action': 'store_true',
                    'help': 'completely skip variant calling, filtering, and annotation',
                }
            },
            {
                'keys': ['--call-region-bed'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'variant calling region bed file (vardict, muse, lofreq, mutect2, haplotype-caller) (default: %(default)s)',
                }
            },
            {
                'keys': ['--panel-of-normal-vcf'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'panel of normal VCF file for Mutect2 (default: %(default)s)',
                }
            },
            {
                'keys': ['--germline-resource-vcf'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'germline resource VCF file for Mutect2 to estimate prior probability (default: %(default)s)',
                }
            },
        ],

    'variant filtering':
        [
            {
                'keys': ['--variant-flagging-criteria'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'comma-separated flagging criteria, e.g. "low_depth: DP<20, mid_qual: 20<=MQ<=40, high_af: AF>0.02" (default: %(default)s)',
                }
            },
            {
                'keys': ['--variant-removal-flags'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'comma-separated flags for variant removal, e.g. "panel_of_normals,map_qual" (default: %(default)s)',
                }
            },
            {
                'keys': ['--only-pass'],
                'properties': {
                    'action': 'store_true',
                    'help': 'only keep the variants with PASS in FILTER column',
                }
            },
        ],

    'variant picking':
        [
            {
                'keys': ['--min-snv-callers'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 1,
                    'help': 'min number of variant callers for an SNV to be picked (default: %(default)s)',
                }
            },
            {
                'keys': ['--min-indel-callers'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 1,
                    'help': 'min number of variant callers for an indel to be picked (default: %(default)s)',
                }
            },
        ],

    'variant annotation':
        [
            {
                'keys': ['--skip-variant-annotation'],
                'properties': {
                    'action': 'store_true',
                    'help': 'do not annotate variants',
                }
            },
            {
                'keys': ['--vep-db-tar-gz'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'VEP database tar.gz file or "vep_cache" dir, required for VEP annotation (default: %(default)s)',
                }
            },
            {
                'keys': ['--vep-db-type'],
                'properties': {
                    'type': str,
                    'required': False,
                    'choices': ['merged', 'vep', 'refseq'],
                    'default': 'merged',
                    'help': 'VEP database type, must match the content in --vep-db-tar-gz (default: %(default)s)',
                }
            },
            {
                'keys': ['--vep-buffer-size'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 5000,
                    'help': 'number of variants loaded in memory by VEP, lower this to reduce memory load (default: %(default)s)',
                }
            },
            {
                'keys': ['--dbnsfp-resource'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'Pre-processed dbNSFP resource file (Bgzip TSV) for VEP (default: %(default)s)',
                }
            },
            {
                'keys': ['--cadd-resource'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'CADD resource file (Bgzip TSV) for VEP (default: %(default)s)',
                }
            },
            {
                'keys': ['--clinvar-vcf-gz'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'ClinVar VCF file (Bgzip VCF) (default: %(default)s)',
                }
            },
            {
                'keys': ['--dbsnp-vcf-gz'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'dbSNP VCF file (Bgzip VCF) (default: %(default)s)',
                }
            },
        ],

    'PCGR':
        [
            {
                'keys': ['--skip-pcgr'],
                'properties': {
                    'action': 'store_true',
                    'help': 'do not run PCGR',
                }
            },
            {
                'keys': ['--pcgr-ref-data-tgz'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'PCGR reference bundle .tgz file (default: %(default)s)',
                }
            },
            {
                'keys': ['--pcgr-vep-tar-gz'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'VEP database tar.gz file or "vep_cache" dir for PCGR (default: %(default)s)',
                }
            },
            {
                'keys': ['--pcgr-tumor-site'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 12,
                    'help': '12 for Head and Neck (default: %(default)s)',
                }
            },
            {
                'keys': ['--pcgr-tmb-target-size-mb'],
                'properties': {
                    'type': int,
                    'required': False,
                    'default': 34,
                    'help': 'effective target size in Mb for TMB analysis (default: %(default)s)',
                }
            },
            {
                'keys': ['--pcgr-tmb-display'],
                'properties': {
                    'type': str,
                    'required': False,
                    'choices': ['coding_and_silent', 'coding_non_silent', 'missense_only'],
                    'default': 'coding_and_silent',
                    'help': 'type of TMB measure to show in report (default: %(default)s)',
                }
            },
        ],

    'vcf2csv':
        [
            {
                'keys': ['-f', '--vcf'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the vcf(.gz) file',
                }
            },
            {
                'keys': ['-c', '--csv'],
                'properties': {
                    'type': str,
                    'required': True,
                    'help': 'path to the output csv file',
                }
            },
        ],

    'others':
        [
            DEBUG_ARG,
            HELP_ARG,
            VERSION_ARG,
        ]
}


class EntryPoint:

    root_parser: argparse.ArgumentParser
    main_parser: argparse.ArgumentParser
    annotate_parser: argparse.ArgumentParser
    vcf2csv_parser: argparse.ArgumentParser

    def main(self):
        self.set_parsers()
        self.add_root_parser_args()
        self.add_main_parser_args()
        self.add_annotate_parser_args()
        self.add_vc2csv_parser_args()
        self.run()

    def set_parsers(self):
        self.root_parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            formatter_class=argparse.RawTextHelpFormatter,
            add_help=False)

        subparsers = self.root_parser.add_subparsers(
            title='commands',
            dest='mode'
        )

        self.main_parser = subparsers.add_parser(
            prog=f'{PROG} main',
            name='main',
            description=f'{DESCRIPTION} - {BOLD}{CYAN}main mode{END}',
            add_help=False)

        self.annotate_parser = subparsers.add_parser(
            prog=f'{PROG} annotate',
            name='annotate',
            description=f'{DESCRIPTION} - {BOLD}{CYAN}annotate mode{END}',
            add_help=False)

        self.vcf2csv_parser = subparsers.add_parser(
            prog=f'{PROG} vcf2csv',
            name='vcf2csv',
            description=f'{DESCRIPTION} - {BOLD}{CYAN}vcf2csv mode{END}',
            add_help=False)

    def add_root_parser_args(self):
        for arg in [HELP_ARG, VERSION_ARG]:
            self.root_parser.add_argument(*arg['keys'], **arg['properties'])

    def add_main_parser_args(self):
        self.__add_arguments(
            parser=self.main_parser,
            required_group_name='main',
            optional_group_names=[
                'normal',
                'general',
                'pre-processing',
                'variant calling',
                'variant filtering',
                'variant picking',
                'variant annotation',
                'PCGR',
            ]
        )

    def add_annotate_parser_args(self):
        self.__add_arguments(
            parser=self.annotate_parser,
            required_group_name='annotate',
            optional_group_names=[
                'general',
                'variant annotation'
            ]
        )

    def add_vc2csv_parser_args(self):
        self.__add_arguments(
            parser=self.vcf2csv_parser,
            required_group_name='vcf2csv',
            optional_group_names=[
                'others',
            ]
        )

    def __add_arguments(
            self,
            parser: argparse.ArgumentParser,
            required_group_name: str,
            optional_group_names: List[str]):

        group = parser.add_argument_group(f'{BOLD}{RED}Required{END}')
        for arg in GROUP_NAME_TO_ARGUMENTS[required_group_name]:
            group.add_argument(*arg['keys'], **arg['properties'])

        for group_name in optional_group_names:
            group = parser.add_argument_group(f'{BOLD}{YELLOW}Optional{END} {BOLD}| {BLUE}{group_name}{END}')
            for arg in GROUP_NAME_TO_ARGUMENTS[group_name]:
                group.add_argument(*arg['keys'], **arg['properties'])

    def run(self):
        args = self.root_parser.parse_args()

        prefix = f'Start running Somatic Pipeline version {__VERSION__} '

        if args.mode is None:
            self.root_parser.print_help()

        elif args.mode == 'main':
            print(f'{prefix}main mode\n', flush=True)
            Run().main(
                ref_fa=args.ref_fa,
                tumor_fq1=args.tumor_fq1,
                tumor_fq2=args.tumor_fq2,

                normal_fq1=args.normal_fq1,
                normal_fq2=args.normal_fq2,
                outdir=args.outdir,
                threads=args.threads,
                debug=args.debug,

                umi_length=args.umi_length,
                clip_r1_5_prime=args.clip_r1_5_prime,
                clip_r2_5_prime=args.clip_r2_5_prime,
                read_aligner=args.read_aligner,
                skip_mark_duplicates=args.skip_mark_duplicates,
                bqsr_known_variant_vcf=args.bqsr_known_variant_vcf,
                discard_bam=args.discard_bam,

                variant_callers=args.variant_callers,
                skip_variant_calling=args.skip_variant_calling,
                call_region_bed=args.call_region_bed,
                panel_of_normal_vcf=args.panel_of_normal_vcf,
                germline_resource_vcf=args.germline_resource_vcf,

                variant_flagging_criteria=args.variant_flagging_criteria,
                variant_removal_flags=args.variant_removal_flags,
                only_pass=args.only_pass,

                min_snv_callers=args.min_snv_callers,
                min_indel_callers=args.min_indel_callers,

                skip_variant_annotation=args.skip_variant_annotation,
                vep_db_tar_gz=args.vep_db_tar_gz,
                vep_db_type=args.vep_db_type,
                vep_buffer_size=args.vep_buffer_size,
                dbnsfp_resource=args.dbnsfp_resource,
                cadd_resource=args.cadd_resource,
                clinvar_vcf_gz=args.clinvar_vcf_gz,
                dbsnp_vcf_gz=args.dbsnp_vcf_gz,

                skip_pcgr=args.skip_pcgr,
                pcgr_ref_data_tgz=args.pcgr_ref_data_tgz,
                pcgr_vep_tar_gz=args.pcgr_vep_tar_gz,
                pcgr_tumor_site=args.pcgr_tumor_site,
                pcgr_tmb_target_size_mb=args.pcgr_tmb_target_size_mb,
                pcgr_tmb_display=args.pcgr_tmb_display
            )

        elif args.mode == 'annotate':
            print(f'{prefix}annotate mode\n', flush=True)
            Run().annotate(
                ref_fa=args.ref_fa,
                vcf=args.vcf,

                outdir=args.outdir,
                threads=args.threads,
                debug=args.debug,

                vep_db_tar_gz=args.vep_db_tar_gz,
                vep_db_type=args.vep_db_type,
                vep_buffer_size=args.vep_buffer_size,
                dbnsfp_resource=args.dbnsfp_resource,
                cadd_resource=args.cadd_resource,
                clinvar_vcf_gz=args.clinvar_vcf_gz,
                dbsnp_vcf_gz=args.dbsnp_vcf_gz
            )

        elif args.mode == 'vcf2csv':
            print(f'{prefix}vcf2csv mode\n', flush=True)
            Run().vcf2csv(
                vcf=args.vcf,
                csv=args.csv,
                debug=args.debug
            )


if __name__ == '__main__':
    EntryPoint().main()
