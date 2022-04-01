import os
from typing import Optional
from .tools import edit_fpath
from .template import Processor


class MappingStats(Processor):

    tumor_bam: str
    normal_bam: Optional[str]

    def main(
            self,
            tumor_bam: str,
            normal_bam: Optional[str]):

        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        SamtoolsStats(self.settings).main(self.tumor_bam)
        if self.normal_bam is not None:
            SamtoolsStats(self.settings).main(self.normal_bam)


class SamtoolsStats(Processor):

    DSTDIR_NAME = 'mapping-stats'

    bam: str

    dstdir: str

    def main(self, bam: str):
        self.bam = bam
        self.make_dstdir()
        self.execute()

    def make_dstdir(self):
        self.dstdir = f'{self.outdir}/{self.DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

    def execute(self):
        txt = edit_fpath(
            fpath=self.bam,
            old_suffix='.bam',
            new_suffix='.txt',
            dstdir=self.dstdir)
        self.call(f'samtools stats {self.bam} > {txt}')
