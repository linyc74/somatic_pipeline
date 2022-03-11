import os
from typing import List, Optional
from .template import Processor


class CNVkit(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    exome_target_bed: str
    gene_annotation_gff: str

    dstdir: str
    mode_args: List[str]

    def main(
            self,
            ref_fa: str,
            tumor_bam: str,
            normal_bam: str,
            gene_annotation_gff: str,
            exome_target_bed: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.gene_annotation_gff = gene_annotation_gff
        self.exome_target_bed = exome_target_bed

        self.make_dstdir()
        self.set_mode_args()
        self.execute()

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/cnvkit'
        os.makedirs(self.dstdir, exist_ok=True)

    def set_mode_args(self):
        if self.exome_target_bed is None:
            mode = 'wgs'
            self.mode_args = [f'-m {mode}']
        else:
            mode = 'hybrid'
            self.mode_args = [
                f'-m {mode}',
                f'--targets {self.exome_target_bed}'
            ]

    def execute(self):
        log = f'{self.outdir}/cnvkit-batch.log'
        args = [
            'cnvkit.py batch',
            self.tumor_bam,
            f'--normal {self.normal_bam}',
            f'--fasta {self.ref_fa}',
            f'--annotate {self.gene_annotation_gff}'
        ] + self.mode_args + [
            '--segment-method cbs',
            f'--output-dir {self.dstdir}',
            f'--processes {self.threads}',
            '--scatter',
            '--diagram',
            f'1> {log} 2> {log}',
        ]
        cmd = self.CMD_LINEBREAK.join(args)
        self.call(cmd)
