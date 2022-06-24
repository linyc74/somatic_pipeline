import argparse
import somatic_pipeline


__VERSION__ = '1.5.2'


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
OPTIONAL_GROUPS = {
    'Optional (general)':
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

    'Optional (pre-processing)':
    [
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

    'Optional (variant calling)':
    [
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
            'keys': ['--skip-variant-calling'],
            'properties': {
                'action': 'store_true',
                'help': 'do not perform variant calling',
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
    ],

    'Optional (variant annotation)':
    [
        {
            'keys': ['--annotator'],
            'properties': {
                'type': str,
                'required': False,
                'choices': ['snpeff', 'vep'],
                'default': 'vep',
                'help': 'variant annotator (default: %(default)s)',
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
                'help': 'ClinVar VCF file (Bgzip VCF), if "None" then skip ClinVar annotation (default: %(default)s)',
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

    'Optional (copy number variation)':
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
        group = self.parser.add_argument_group('Required')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        for group_name, arguments in OPTIONAL_GROUPS.items():
            group = self.parser.add_argument_group(group_name)
            for arg in arguments:
                group.add_argument(*arg['keys'], **arg['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running Somatic Pipeline version {__VERSION__}\n', flush=True)
        somatic_pipeline.Main().main(
            ref_fa=args.ref_fa,
            tumor_fq1=args.tumor_fq1,
            tumor_fq2=args.tumor_fq2,

            normal_fq1=args.normal_fq1,
            normal_fq2=args.normal_fq2,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug,

            read_aligner=args.read_aligner,
            skip_mark_duplicates=args.skip_mark_duplicates,
            bqsr_known_variant_vcf=args.bqsr_known_variant_vcf,
            discard_bam=args.discard_bam,

            variant_caller=args.variant_caller,
            skip_variant_calling=args.skip_variant_calling,
            panel_of_normal_vcf=args.panel_of_normal_vcf,

            annotator=args.annotator,
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
            cnvkit_annotate_txt=args.cnvkit_annotate_txt)


if __name__ == '__main__':
    EntryPoint().main()
