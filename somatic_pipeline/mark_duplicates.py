from typing import Optional, Tuple
from .tools import edit_fpath
from .template import Processor


class MarkDuplicates(Processor):

    tumor_bam: str
    normal_bam: Optional[str]

    out_tumor_bam: str
    out_normal_bam: Optional[str]

    def main(
            self,
            tumor_bam: str,
            normal_bam: Optional[str]) -> Tuple[str, Optional[str]]:

        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        self.out_tumor_bam = GATKMarkDuplicates(self.settings).main(
            bam=self.tumor_bam)

        self.out_normal_bam = None if self.normal_bam is None else \
            GATKMarkDuplicates(self.settings).main(bam=self.normal_bam)

        return self.out_tumor_bam, self.out_normal_bam


class GATKMarkDuplicates(Processor):

    REMOVE_DUPLICATES = 'true'

    bam: str

    metrics_txt: str
    out_bam: str

    def main(self, bam: str) -> str:
        self.bam = bam
        self.set_out_bam()
        self.set_metrics_txt()
        self.execute()
        return self.out_bam

    def set_out_bam(self):
        self.out_bam = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='-mark-duplicates.bam',
            dstdir=self.workdir)

    def set_metrics_txt(self):
        self.metrics_txt = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='-duplicate-metrics.txt',
            dstdir=self.outdir)

    def execute(self):
        log = f'{self.outdir}/gatk-MarkDuplicates.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk MarkDuplicates',
            f'--INPUT {self.bam}',
            f'--METRICS_FILE {self.metrics_txt}',
            f'--OUTPUT {self.out_bam}',
            f'--REMOVE_DUPLICATES {self.REMOVE_DUPLICATES}',
            f'1> {log} 2> {log}',
        ])
        self.call(cmd)
