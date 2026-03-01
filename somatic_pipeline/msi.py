import os
from typing import Optional
from .template import Processor
from .index_files import SamtoolsIndexBam


class MSIsensor(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: str
    bed_file: str

    def main(self, ref_fa: str, tumor_bam: str, normal_bam: str, bed_file: Optional[str]):

        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.bed_file = bed_file

        lines = [
            'msisensor scan',
            f'-d {self.ref_fa}',
            f'-o {self.workdir}/microsatellites.list',
            f'1> {self.outdir}/msisensor-scan.log',
            f'2> {self.outdir}/msisensor-scan.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        for bam in [self.tumor_bam, self.normal_bam]:
            if not os.path.exists(f'{bam}.bai'):
                SamtoolsIndexBam(self.settings).main(bam=bam)

        os.makedirs(f'{self.outdir}/msi', exist_ok=True)

        lines = [
            ' msisensor msi',
            f'-d {self.workdir}/microsatellites.list',
            f'-n {self.normal_bam}',
            f'-t {self.tumor_bam}',
        ]
        if self.bed_file is not None:
            lines += [f'-e {self.bed_file}']
        lines += [
            f'-o {self.outdir}/msi/msisensor',
            f'1> {self.outdir}/msisensor-msi.log',
            f'2> {self.outdir}/msisensor-msi.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))


class MANTIS(Processor):

    ref_fa: str
    tumor_bam: str
    normal_bam: str

    def main(self, ref_fa: str, tumor_bam: str, normal_bam: str):
        self.ref_fa = ref_fa
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam

        lines = [
            f'RepeatFinder',
            f'-i {self.ref_fa}',
            f'-o {self.workdir}/RepeatFinder-microsatellites.bed',
            f'1> {self.outdir}/RepeatFinder.log',
            f'2> {self.outdir}/RepeatFinder.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))

        for bam in [self.tumor_bam, self.normal_bam]:
            if not os.path.exists(f'{bam}.bai'):
                SamtoolsIndexBam(self.settings).main(bam=bam)
        
        os.makedirs(f'{self.outdir}/msi', exist_ok=True)

        lines = [
            f'mantis.py',
            f'--bedfile {self.workdir}/RepeatFinder-microsatellites.bed',
            f'--genome {self.ref_fa}',
            f'-n {self.normal_bam}',
            f'-t {self.tumor_bam}',
            f'--threads {self.threads}',
            f'--min-read-quality 20.0',
            f'--min-locus-quality 25.0',
            f'--min-locus-coverage 20',
            f'--min-repeat-reads 1',
            f'-o {self.outdir}/msi/mantis.txt',
            f'1> {self.outdir}/mantis.log',
            f'2> {self.outdir}/mantis.log',
        ]
        self.call(self.CMD_LINEBREAK.join(lines))
