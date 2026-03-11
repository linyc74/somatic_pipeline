import os
from os.path import basename
from typing import List, Optional
from .tools import edit_fpath
from .template import Processor


class CNV(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    exome_target_bed: Optional[str]
    annotate_txt: Optional[str]
    segmentation_threshold: float

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            exome_target_bed: Optional[str],
            annotate_txt: Optional[str],
            segmentation_threshold: float):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.exome_target_bed = exome_target_bed
        self.annotate_txt = annotate_txt
        self.segmentation_threshold = segmentation_threshold
        
        if self.exome_target_bed is not None:
            self.exome_target_bed = CleanUpBed(self.settings).main(
                bed=self.exome_target_bed)

        CNVkitBatch(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            exome_target_bed=self.exome_target_bed,
            annotate_txt=self.annotate_txt,
            segmentation_threshold=self.segmentation_threshold)


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


class CNVkitBatch(Processor):

    SEGMENT_METHOD = 'cbs'  # depends on the R package DNAcopy

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    exome_target_bed: Optional[str]
    annotate_txt: Optional[str]
    segmentation_threshold: float

    dstdir: str
    cnr_file: str  # bin-level copy number ratio
    cns_file: str  # copy number segmentation

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            exome_target_bed: Optional[str],
            annotate_txt: Optional[str],
            segmentation_threshold: float):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.exome_target_bed = exome_target_bed
        self.annotate_txt = annotate_txt
        self.segmentation_threshold = segmentation_threshold
        
        self.dstdir = f'{self.workdir}/cnvkit'
        os.makedirs(self.dstdir, exist_ok=True)

        self.cnvkit_batch()
        self.set_cnr_cns_filepaths()
        self.cnvkit_segment()
        self.cnvkit_plots()
        self.move_files_to_outdir()

    def cnvkit_batch(self):
        lines = [
            'cnvkit.py batch',
            f'--normal {self.normal_bam}',
            f'--fasta {self.ref_fa}',
            f'--segment-method {self.SEGMENT_METHOD}',
            '--drop-low-coverage',
            f'--output-dir {self.dstdir}',
            f'--processes {self.threads}',
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

    def set_cnr_cns_filepaths(self):
        fname = basename(self.tumor_bam)[:-len('.bam')]
        self.cnr_file = f'{self.dstdir}/{fname}.cnr'
        self.cns_file = f'{self.dstdir}/{fname}.cns'

    def cnvkit_segment(self):
        # the `batch` command already runs the `segment` command internally
        # this is to rerun the `segment` command with a different threshold
        # to override the original segmentation results
        lines = [
            'cnvkit.py segment',
            self.cnr_file,
            f'--method {self.SEGMENT_METHOD}',
            f'--threshold {self.segmentation_threshold}',
            '--drop-low-coverage',
            f'--output {self.cns_file}',
            f'1> {self.outdir}/cnvkit-segment.log',
            f'2> {self.outdir}/cnvkit-segment.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
    
    def cnvkit_plots(self):
        lines = [
            'cnvkit.py scatter',
            self.cnr_file,
            f'--segment {self.cns_file}',
            f'--output {self.dstdir}/scatter.pdf',
            f'1> {self.outdir}/cnvkit-scatter.log',
            f'2> {self.outdir}/cnvkit-scatter.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        lines = [
            'cnvkit.py diagram',
            self.cns_file,
            f'--output {self.dstdir}/diagram.pdf',
            f'1> {self.outdir}/cnvkit-diagram.log',
            f'2> {self.outdir}/cnvkit-diagram.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
    
    def move_files_to_outdir(self):
        os.makedirs(f'{self.outdir}/cnvkit', exist_ok=True)
        for f in [
            self.cnr_file,
            self.cns_file,
            f'{self.dstdir}/scatter.pdf',
            f'{self.dstdir}/diagram.pdf',
        ]:
            self.call(f'mv {f} {self.outdir}/cnvkit/')
