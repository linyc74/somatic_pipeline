import argparse
import gatk_pipeline


__VERSION__ = '1.1.0-beta'


PROG = 'python gatk_pipeline'
DESCRIPTION = f'Custom-built GATK pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-r', '--ref-fa'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the reference genome fasta file',
        }
    },
    {
        'keys': ['-1', '--tumor-fq1'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the tumor read 1 fastq file',
        }
    },
    {
        'keys': ['-2', '--tumor-fq2'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the tumor read 2 fastq file',
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
            'help': 'path to the normal read 1 fastq file (default: %(default)s)',
        }
    },
    {
        'keys': ['-4', '--normal-fq2'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'path to the normal read 2 fastq file (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'gatk_pipeline_outdir',
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
        print(f'Start running GATK Pipeline version {__VERSION__}\n', flush=True)
        gatk_pipeline.Main().main(
            ref_fa=args.ref_fa,
            tumor_fq1=args.tumor_fq1,
            tumor_fq2=args.tumor_fq2,
            normal_fq1=args.normal_fq1,
            normal_fq2=args.normal_fq2,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
