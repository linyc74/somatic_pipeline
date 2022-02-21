import os
from typing import Optional
from .template import Settings
from .somatic_pipeline import SomaticPipeline


class Main:

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    read_aligner: str
    variant_caller: str

    settings: Settings

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: str,
            normal_fq2: str,
            read_aligner: str,
            variant_caller: str,
            outdir: str,
            threads: str,
            debug: bool):

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = None if normal_fq1 == 'None' else normal_fq1
        self.normal_fq2 = None if normal_fq2 == 'None' else normal_fq2
        self.read_aligner = read_aligner
        self.variant_caller = variant_caller

        self.settings = Settings(
            workdir='./gatk_pipeline_workdir',
            outdir=outdir,
            threads=int(threads),
            debug=debug,
            mock=False)

        for d in [self.settings.workdir, self.settings.outdir]:
            os.makedirs(d, exist_ok=True)

        SomaticPipeline(self.settings).main(
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,
            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2,
            read_aligner=self.read_aligner,
            variant_caller=self.variant_caller)
