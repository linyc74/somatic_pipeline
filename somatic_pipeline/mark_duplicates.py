from .tools import edit_fpath
from .template import Processor


class MarkDuplicates(Processor):

    bam_in: str

    metrics_txt: str
    bam_out: str

    def main(self, bam_in: str) -> str:
        self.bam_in = bam_in
        self.set_bam_out()
        self.set_metrics_txt()
        self.execute()
        return self.bam_out

    def set_bam_out(self):
        self.bam_out = edit_fpath(
            fpath=self.bam_in,
            old_suffix='.bam',
            new_suffix='-mark-duplicates.bam',
            dstdir=self.outdir)

    def set_metrics_txt(self):
        self.metrics_txt = edit_fpath(
            fpath=self.bam_in,
            old_suffix='.bam',
            new_suffix='-duplicate-metrics.txt',
            dstdir=self.outdir)

    def execute(self):
        log = f'{self.outdir}/gatk-MarkDuplicates.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk MarkDuplicates',
            f'--INPUT {self.bam_in}',
            f'--METRICS_FILE {self.metrics_txt}',
            f'--OUTPUT {self.bam_out}',
            f'1> {log} 2> {log}',
        ])
        self.call(cmd)
