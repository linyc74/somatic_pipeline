from typing import Optional, Tuple
from abc import ABC, abstractmethod
from os.path import exists, dirname, samefile
from .template import Processor
from .constant import TUMOR, NORMAL


class Mapping(Processor):

    read_aligner: str
    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    discard_bam: bool

    index: str
    tumor_bam: str
    normal_bam: Optional[str]

    def main(
            self,
            read_aligner: str,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str],
            discard_bam: bool) -> Tuple[str, Optional[str]]:

        self.read_aligner = read_aligner
        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.discard_bam = discard_bam

        if self.read_aligner == 'bwa':
            indexer = BwaIndexer(self.settings).main
            aligner = BwaAligner(self.settings).main
        elif self.read_aligner == 'bowtie2':
            indexer = Bowtie2Indexer(self.settings).main
            aligner = Bowtie2Aligner(self.settings).main
        else:
            raise ValueError(f'Invalid read aligner: {self.read_aligner}')

        self.index = indexer(fna=self.ref_fa)

        self.tumor_bam = aligner(
            index=self.index,
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2,
            sample_name=TUMOR,
            discard_bam=self.discard_bam)

        if self.normal_fq1 is not None:
            self.normal_bam = aligner(
                index=self.index,
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                sample_name=NORMAL,
                discard_bam=self.discard_bam)
        else:
            self.normal_bam = None

        self.call(f'rm {self.index}*')  # to save disk space

        return self.tumor_bam, self.normal_bam


class Indexer(Processor, ABC):

    fna: str

    index: str

    def main(self, fna: str) -> str:
        self.fna = fna
        self.run_indexing()
        return self.index

    @abstractmethod
    def run_indexing(self):
        return


class BwaIndexer(Indexer):

    def run_indexing(self):

        self.index = f'{self.workdir}/bwa-index'

        indexed = True
        for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
            if not exists(f'{self.index}{ext}'):
                indexed = False
        if indexed:
            return

        log = f'{self.outdir}/bwa-index.log'
        cmd = self.CMD_LINEBREAK.join([
            'bwa index',
            f'-p {self.index}',
            self.fna,
            f'2> {log}',
        ])
        self.call(cmd)


class Bowtie2Indexer(Indexer):

    def run_indexing(self):
        self.index = f'{self.workdir}/bowtie2-index'
        log = f'{self.outdir}/bowtie2-build.log'
        cmd = self.CMD_LINEBREAK.join([
            'bowtie2-build',
            self.fna,
            self.index,
            f'1> {log}',
            f'2> {log}',
        ])
        self.call(cmd)


class Aligner(Processor, ABC):

    index: str
    fq1: str
    fq2: str
    sample_name: str
    discard_bam: bool

    sam: str
    bam: str
    sorted_bam: str

    def main(
            self,
            index: str,
            fq1: str,
            fq2: str,
            sample_name: str,
            discard_bam: bool) -> str:

        self.index = index
        self.fq1 = fq1
        self.fq2 = fq2
        self.sample_name = sample_name
        self.discard_bam = discard_bam

        self.align()
        self.remove_fastqs()
        self.sam_to_bam()
        self.sort_bam()

        return self.sorted_bam

    @abstractmethod
    def align(self):
        return

    def remove_fastqs(self):
        for fq in [self.fq1, self.fq2]:
            if samefile(dirname(fq), self.workdir):  # fq in the workdir
                self.call(f'rm {fq}')

    def sam_to_bam(self):
        self.bam = f'{self.workdir}/{self.sample_name}.bam'
        self.call(f'samtools view -b -h {self.sam} > {self.bam}')
        self.call(f'rm {self.sam}')  # to save disk space

    def sort_bam(self):
        dstdir = self.workdir if self.discard_bam else self.outdir
        self.sorted_bam = f'{dstdir}/{self.sample_name}-sorted.bam'
        self.call(f'samtools sort {self.bam} > {self.sorted_bam}')
        self.call(f'rm {self.bam}')  # to save disk space


class BwaAligner(Aligner):

    def align(self):
        self.sam = f'{self.workdir}/{self.sample_name}.sam'
        log = f'{self.outdir}/bwa-mem-{self.sample_name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'bwa mem',
            f'-t {self.threads}',
            f'-o {self.sam}',
            f'-R "@RG\\tID:foo\\tSM:{self.sample_name}\\tPL:ILLUMINA\\tLB:{self.sample_name}"',  # needs literal '\t' for tab
            self.index,
            self.fq1,
            self.fq2,
            f'2> {log}',
        ])
        self.call(cmd)


class Bowtie2Aligner(Aligner):

    def align(self):
        self.sam = f'{self.workdir}/{self.sample_name}.sam'
        log = f'{self.outdir}/bowtie2-{self.sample_name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'bowtie2',
            f'-x {self.index}',
            f'-1 {self.fq1}',
            f'-2 {self.fq2}',
            f'-S {self.sam}',
            f'--rg-id foo',
            f'--rg "SM:{self.sample_name}"',
            f'--rg "PL:ILLUMINA"',
            f'--rg "LB:{self.sample_name}"',  # library tag required by some callers, e.g. SomaticSniper
            f'--threads {self.threads}',
            f'2> {log}',
        ])
        self.call(cmd)
