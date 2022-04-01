import argparse
import somatic_pipeline


__VERSION__ = '1.3.4-beta'


PROG = 'python somatic_pipeline'
DESCRIPTION = f'Custom-built somatic pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
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
]
OPTIONAL = [
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
        'keys': ['--variant-caller'],
        'properties': {
            'type': str,
            'required': False,
            'choices': ['mutect2', 'muse', 'varscan'],
            'default': 'mutect2',
            'help': 'variant caller (default: %(default)s)',
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
        'keys': ['--bqsr-known-variant-vcf'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'known variants VCF file for BQSR, if None then skip BQSR (default: %(default)s)',
        }
    },
    {
        'keys': ['--discard-bam'],
        'properties': {
            'action': 'store_true',
            'help': 'do not save sorted BAM files in outdir',
        }
    },
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
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running Somatic Pipeline version {__VERSION__}\n', flush=True)
        somatic_pipeline.Main().main(
            ref_fa=args.ref_fa,
            tumor_fq1=args.tumor_fq1,
            tumor_fq2=args.tumor_fq2,
            normal_fq1=args.normal_fq1,
            normal_fq2=args.normal_fq2,
            read_aligner=args.read_aligner,
            variant_caller=args.variant_caller,
            exome_target_bed=args.exome_target_bed,
            cnvkit_annotate_txt=args.cnvkit_annotate_txt,
            panel_of_normal_vcf=args.panel_of_normal_vcf,
            bqsr_known_variant_vcf=args.bqsr_known_variant_vcf,
            discard_bam=args.discard_bam,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
