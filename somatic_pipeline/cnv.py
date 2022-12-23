import os
from typing import List, Optional
from .tools import edit_fpath
from .template import Processor


class ComputeCNV(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: Optional[str]
    exome_target_bed: Optional[str]
    annotate_txt: Optional[str]

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: Optional[str],
            exome_target_bed: Optional[str],
            annotate_txt: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.exome_target_bed = exome_target_bed
        self.annotate_txt = annotate_txt

        if self.normal_bam is None:
            self.logger.info(f'Normal sample not provided, skip CNV calculation')
            return

        self.clean_up_bed()
        self.report_auto_bin_size()
        self.run_cnvkit()

    def clean_up_bed(self):
        if self.exome_target_bed is not None:
            self.exome_target_bed = CleanUpBed(self.settings).main(
                bed=self.exome_target_bed)

    def report_auto_bin_size(self):
        if self.exome_target_bed is not None:
            CNVkitReportAutoBinSize(self.settings).main(
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                exome_target_bed=self.exome_target_bed)

    def run_cnvkit(self):
        CNVkitBatch(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            exome_target_bed=self.exome_target_bed,
            annotate_txt=self.annotate_txt)


class CleanUpBed(Processor):

    bed: str

    output_bed: str

    def main(self, bed: str) -> str:
        self.bed = bed
        self.set_output_bed()
        self.clean_up()
        return self.output_bed

    def set_output_bed(self):
        self.output_bed = edit_fpath(
            fpath=self.bed,
            old_suffix='.call_region_bed',
            new_suffix='-clean.call_region_bed',
            dstdir=self.workdir)

    def clean_up(self):
        with open(self.bed) as reader:
            with open(self.output_bed, 'w') as writer:
                for line in reader:
                    fields = len(line.strip().split('\t'))
                    if fields >= 3:
                        writer.write(line)


class CNVkitBatch(Processor):

    DSTDIR_NAME = 'cnvkit'
    SEGMENT_METHOD = 'cbs'

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    exome_target_bed: Optional[str]
    annotate_txt: Optional[str]

    dstdir: str
    method_args: List[str]
    annotate_args: List[str]

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            exome_target_bed: Optional[str],
            annotate_txt: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.exome_target_bed = exome_target_bed
        self.annotate_txt = annotate_txt

        self.make_dstdir()
        self.set_method_args()
        self.set_annotate_args()
        self.execute()

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def set_method_args(self):
        if self.exome_target_bed is None:
            self.method_args = [f'--method wgs']
        else:
            self.method_args = [
                f'--method hybrid',
                f'--targets {self.exome_target_bed}'
            ]

    def set_annotate_args(self):
        self.annotate_args = [] \
            if self.annotate_txt is None \
            else [f'--annotate {self.annotate_txt}']

    def execute(self):
        log = f'{self.outdir}/cnvkit-batch.log'
        args = [
            'cnvkit.py batch',
            f'--normal {self.normal_bam}',
            f'--fasta {self.ref_fa}',
        ] + self.annotate_args + self.method_args + [
            f'--segment-method {self.SEGMENT_METHOD}',
            '--drop-low-coverage',
            f'--output-dir {self.dstdir}',
            f'--processes {self.threads}',
            '--scatter',
            '--diagram',
            self.tumor_bam,
            f'1> {log}',
            f'2> {log}',
        ]
        cmd = self.CMD_LINEBREAK.join(args)
        self.call(cmd)


class CNVkitReportAutoBinSize(Processor):

    tumor_bam: str
    normal_bam: str
    exome_target_bed: str

    def main(
            self,
            tumor_bam: str,
            normal_bam: str,
            exome_target_bed: str):

        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.exome_target_bed = exome_target_bed

        self.execute()

    def execute(self):
        log = f'{self.outdir}/cnvkit-autobin.log'
        cmd = self.CMD_LINEBREAK.join([
            'cnvkit.py autobin',
            self.normal_bam,
            self.tumor_bam,
            '--method hybrid',
            f'--targets {self.exome_target_bed}',
            f'--target-output-call_region_bed {self.workdir}/cnvkit-autobin-target.call_region_bed',
            f'--antitarget-output-call_region_bed {self.workdir}/cnvkit-autobin-antitarget.call_region_bed',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)
