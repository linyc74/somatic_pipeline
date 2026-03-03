import os
from typing import List, Optional
from .tools import edit_fpath
from .template import Processor


class CNV(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    exome_target_bed: Optional[str]
    annotate_txt: Optional[str]

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

        if self.exome_target_bed is not None:
            self.exome_target_bed = CleanUpBed(self.settings).main(
                bed=self.exome_target_bed)

        if self.exome_target_bed is not None:
            CNVkitReportAutoBinSize(self.settings).main(
                tumor_bam=self.tumor_bam,
                normal_bam=self.normal_bam,
                exome_target_bed=self.exome_target_bed)

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
            old_suffix='.bed',
            new_suffix='-clean.bed',
            dstdir=self.workdir)

    def clean_up(self):
        with open(self.bed) as reader:
            with open(self.output_bed, 'w') as writer:
                for line in reader:
                    fields = len(line.strip().split('\t'))
                    if fields >= 3:
                        writer.write(line)


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
            f'--target-output-bed {self.workdir}/cnvkit-autobin-target.bed',
            f'--antitarget-output-bed {self.workdir}/cnvkit-autobin-antitarget.bed',
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class CNVkitBatch(Processor):

    DSTDIR_NAME = 'cnvkit'
    SEGMENT_METHOD = 'cbs'  # depends on the R package DNAcopy

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    exome_target_bed: Optional[str]
    annotate_txt: Optional[str]

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

        dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(dstdir, exist_ok=True)

        lines = [
            'cnvkit.py batch',
            f'--normal {self.normal_bam}',
            f'--fasta {self.ref_fa}',
            f'--segment-method {self.SEGMENT_METHOD}',
            '--drop-low-coverage',
            f'--output-dir {dstdir}',
            f'--processes {self.threads}',
            '--scatter',
            '--diagram',
        ]

        if self.exome_target_bed is not None:
            lines += [
                f'--method hybrid',  # hybrid (WES) mode requires target BED file
                f'--targets {self.exome_target_bed}'
            ]
        else:  # if no target BED file, can only run WGS mode
            lines += [f'--method wgs']
        
        if self.annotate_txt is not None:
            lines += [f'--annotate {self.annotate_txt}']

        lines += [
            self.tumor_bam,
            f'1> {self.outdir}/cnvkit-batch.log',
            f'2> {self.outdir}/cnvkit-batch.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
