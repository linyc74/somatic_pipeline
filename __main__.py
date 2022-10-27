import argparse
from somatic_pipeline import Run


__VERSION__ = '1.7.0-beta'


PROG = 'python somatic_pipeline'
DESCRIPTION = f'Somatic pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED_GROUP_NAME_TO_ARGUMENTS = {
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
}
OPTIONAL_GROUP_NAME_TO_ARGUMENTS = {
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
            {
                'keys': ['-d', '--debug'],
                'properties': {
                    'action': 'store_true',
                    'help': 'debug mode',
                }
            },
            {
                'keys': ['-h', '--help'],
                'properties': {
                    'action': 'help',
                    'help': 'show this help message',
                }
            },
            {
                'keys': ['-v', '--version'],
                'properties': {
                    'action': 'version',
                    'version': __VERSION__,
                    'help': 'show version',
                }
            },
        ],

    'pre-processing':
        [
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
                'keys': ['--variant-caller'],
                'properties': {
                    'type': str,
                    'required': False,
                    'choices': ['mutect2', 'muse', 'varscan', 'haplotype-caller', 'vardict', 'lofreq', 'somatic-sniper'],
                    'default': 'mutect2',
                    'help': '"mutect2", "haplotype-caller", "vardict", "lofreq" works for tumor-only mode (default: %(default)s)',
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
            {
                'keys': ['--vardict-call-region-bed'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'variant calling region bed file (WES) for VarDict (default: %(default)s)',
                }
            },
        ],

    'variant filtering':
        [
            {
                'keys': ['--variant-removal-flags'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'comma-separated flags for variant removal, e.g. "panel_of_normals,map_qual" (default: %(default)s)',
                }
            },
        ],

    'variant annotation':
        [
            {
                'keys': ['--annotator'],
                'properties': {
                    'type': str,
                    'required': False,
                    'choices': ['vep', 'snpeff'],
                    'default': 'vep',
                    'help': 'variant annotator (default: %(default)s)',
                }
            },
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
                    'help': 'VEP database tar.gz file, required for VEP annotation (default: %(default)s)',
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
            {
                'keys': ['--snpsift-dbnsfp-txt-gz'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'SnpSift dbNSFP database file (Bgzip VCF) (default: %(default)s)',
                }
            },
        ],

    'copy number variation':
        [
            {
                'keys': ['--skip-cnv'],
                'properties': {
                    'action': 'store_true',
                    'help': 'do not compute copy number variation',
                }
            },
            {
                'keys': ['--exome-target-bed'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'BED file of exome target probes to compute CNV (default: %(default)s)',
                }
            },
            {
                'keys': ['--cnvkit-annotate-txt'],
                'properties': {
                    'type': str,
                    'required': False,
                    'default': 'None',
                    'help': 'CNVkit annotation file, usually UCSC refFlat.txt (default: %(default)s)',
                }
            },
        ],
}


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


class EntryPoint:

    __root: argparse.ArgumentParser
    __main: argparse.ArgumentParser
    __annotate: argparse.ArgumentParser

    def main(self):
        self.set_parsers()
        self.add_main_parser_arguments()
        self.add_annotate_parser_arguments()
        self.run()

    def set_parsers(self):
        self.__root = argparse.ArgumentParser(
            description=DESCRIPTION,
            formatter_class=argparse.RawTextHelpFormatter)

        subparsers = self.__root.add_subparsers(
            title='commands',
            dest='mode'
        )

        self.__main = subparsers.add_parser(
            prog=PROG,
            name='main',
            description=DESCRIPTION,
            add_help=False)

        self.__annotate = subparsers.add_parser(
            prog=PROG,
            name='annotate',
            description=DESCRIPTION,
            add_help=False)

    def add_main_parser_arguments(self):
        group = self.__main.add_argument_group(f'{BOLD}{RED}Required{END}')
        for item in REQUIRED_GROUP_NAME_TO_ARGUMENTS['main']:
            group.add_argument(*item['keys'], **item['properties'])

        for group_name, arguments in OPTIONAL_GROUP_NAME_TO_ARGUMENTS.items():
            group = self.__main.add_argument_group(f'{BOLD}{YELLOW}Optional {BLUE}({group_name}){END}')
            for arg in arguments:
                group.add_argument(*arg['keys'], **arg['properties'])

    def add_annotate_parser_arguments(self):
        g = self.__annotate.add_argument_group(f'{BOLD}{RED}Required{END}')
        for item in REQUIRED_GROUP_NAME_TO_ARGUMENTS['annotate']:
            g.add_argument(*item['keys'], **item['properties'])

        for group_name in ['general', 'variant annotation']:
            g = self.__annotate.add_argument_group(f'{BOLD}{YELLOW}Optional {BLUE}({group_name}){END}')
            arguments = OPTIONAL_GROUP_NAME_TO_ARGUMENTS[group_name]
            for arg in arguments:
                g.add_argument(*arg['keys'], **arg['properties'])

    def run(self):
        args = self.__root.parse_args()

        if args.mode == 'main':
            print(f'Start running Somatic Pipeline version {__VERSION__}\n', flush=True)
            Run().main(
                ref_fa=args.ref_fa,
                tumor_fq1=args.tumor_fq1,
                tumor_fq2=args.tumor_fq2,

                normal_fq1=args.normal_fq1,
                normal_fq2=args.normal_fq2,
                outdir=args.outdir,
                threads=args.threads,
                debug=args.debug,

                clip_r1_5_prime=args.clip_r1_5_prime,
                clip_r2_5_prime=args.clip_r2_5_prime,
                read_aligner=args.read_aligner,
                skip_mark_duplicates=args.skip_mark_duplicates,
                bqsr_known_variant_vcf=args.bqsr_known_variant_vcf,
                discard_bam=args.discard_bam,

                variant_caller=args.variant_caller,
                skip_variant_calling=args.skip_variant_calling,
                panel_of_normal_vcf=args.panel_of_normal_vcf,
                germline_resource_vcf=args.germline_resource_vcf,
                vardict_call_region_bed=args.vardict_call_region_bed,

                variant_removal_flags=args.variant_removal_flags,

                annotator=args.annotator,
                skip_variant_annotation=args.skip_variant_annotation,
                vep_db_tar_gz=args.vep_db_tar_gz,
                vep_db_type=args.vep_db_type,
                vep_buffer_size=args.vep_buffer_size,
                dbnsfp_resource=args.dbnsfp_resource,
                cadd_resource=args.cadd_resource,
                clinvar_vcf_gz=args.clinvar_vcf_gz,
                dbsnp_vcf_gz=args.dbsnp_vcf_gz,
                snpsift_dbnsfp_txt_gz=args.snpsift_dbnsfp_txt_gz,

                skip_cnv=args.skip_cnv,
                exome_target_bed=args.exome_target_bed,
                cnvkit_annotate_txt=args.cnvkit_annotate_txt
            )

        elif args.mode == 'annotate':
            print(f'Start running Somatic Pipeline version {__VERSION__} \033[93mannotate\033[0m\n', flush=True)
            Run().annotate(
                ref_fa=args.ref_fa,
                vcf=args.vcf,

                outdir=args.outdir,
                threads=args.threads,
                debug=args.debug,

                annotator=args.annotator,
                vep_db_tar_gz=args.vep_db_tar_gz,
                vep_db_type=args.vep_db_type,
                vep_buffer_size=args.vep_buffer_size,
                dbnsfp_resource=args.dbnsfp_resource,
                cadd_resource=args.cadd_resource,
                clinvar_vcf_gz=args.clinvar_vcf_gz,
                dbsnp_vcf_gz=args.dbsnp_vcf_gz,
                snpsift_dbnsfp_txt_gz=args.snpsift_dbnsfp_txt_gz
            )


if __name__ == '__main__':
    EntryPoint().main()
