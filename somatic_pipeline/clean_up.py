import os
from .template import Processor, Settings


class CleanUp(Processor):

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self):
        self.collect_log_files()
        self.remove_workdir()

    def collect_log_files(self):
        os.makedirs(f'{self.outdir}/log')
        cmd = f'mv {self.outdir}/*.log {self.outdir}/log/'
        self.call(cmd)

    def remove_workdir(self):
        if not self.debug:
            self.call(f'rm -r {self.workdir}')
