from typing import Optional, Tuple
from abc import ABC, abstractmethod
from .constant import TUMOR, NORMAL
from .template import Processor, Settings


class TemplateIndexer(Processor, ABC):

    fna: str
    index: str

    def main(self, fna: str) -> str:
        self.fna = fna
        self.run_indexing()
        return self.index

    @abstractmethod
    def run_indexing(self):
        return


class BwaIndexer(TemplateIndexer):

    def run_indexing(self):
        self.index = f'{self.workdir}/bwa-index'
        log = f'{self.outdir}/bwa-index.log'
        cmd = self.CMD_LINEBREAK.join([
            'bwa index',
            f'-p {self.index}',
            self.fna,
            f'2> {log}',
        ])
        self.call(cmd)


class Bowtie2Indexer(TemplateIndexer):

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


class TemplateAligner(Processor, ABC):

    index: str
    fq1: str
    fq2: str
    sample_name: str

    sam: str
    bam: str
    sorted_bam: str

    def main(
            self,
            index: str,
            fq1: str,
            fq2: str,
            sample_name: str) -> str:

        self.index = index
        self.fq1 = fq1
        self.fq2 = fq2
        self.sample_name = sample_name

        self.align()
        self.sam_to_bam()
        self.sort_bam()

        return self.sorted_bam

    @abstractmethod
    def align(self):
        return

    def sam_to_bam(self):
        self.bam = f'{self.workdir}/{self.sample_name}.bam'
        cmd = f'samtools view -b -h {self.sam} > {self.bam}'
        self.call(cmd)

    def sort_bam(self):
        self.sorted_bam = f'{self.outdir}/{self.sample_name}-sorted.bam'
        cmd = f'samtools sort {self.bam} > {self.sorted_bam}'
        self.call(cmd)


class BwaAligner(TemplateAligner):

    def align(self):
        self.sam = f'{self.workdir}/{self.sample_name}.sam'
        log = f'{self.outdir}/bwa-mem-{self.sample_name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'bwa mem',
            f'-t {self.threads}',
            f'-o {self.sam}',
            self.index,
            self.fq1,
            self.fq2,
            f'2> {log}',
        ])
        self.call(cmd)


class Bowtie2Aligner(TemplateAligner):

    def align(self):
        self.sam = f'{self.workdir}/{self.sample_name}.sam'
        log = f'{self.outdir}/bowtie2-{self.sample_name}.log'
        cmd = self.CMD_LINEBREAK.join([
            'bowtie2',
            f'-x {self.index}',
            f'-1 {self.fq1}',
            f'-2 {self.fq2}',
            f'-S {self.sam}',
            f'--threads {self.threads}',
            f'2> {log}',
        ])
        self.call(cmd)


class StrategyIndexAndAlignTumorNormal:

    indexer: TemplateIndexer
    aligner: TemplateAligner

    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]

    index: str
    tumor_bam: str
    normal_bam: Optional[str]

    def __init__(self, indexer: TemplateIndexer, aligner: TemplateAligner):
        self.indexer = indexer
        self.aligner = aligner

    def main(
            self,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str]) -> Tuple[str, Optional[str]]:

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2

        self.index_genome()
        self.align_tumor()
        self.align_normal()

        return self.tumor_bam, self.normal_bam

    def index_genome(self):
        self.index = self.indexer.main(fna=self.ref_fa)

    def align_tumor(self):
        self.tumor_bam = self.aligner.main(
            index=self.index,
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2,
            sample_name=TUMOR)

    def align_normal(self):
        self.normal_bam = None if self.normal_fq1 is None \
            else self.aligner.main(
                index=self.index,
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                sample_name=NORMAL)


class FactoryIndexAndAlignTumorNormal:

    INDEXER_DICT = {
        'bwa': BwaIndexer,
        'bowtie2': Bowtie2Indexer,
    }
    ALIGNER_DICT = {
        'bwa': BwaAligner,
        'bowtie2': Bowtie2Aligner,
    }

    def get_callable(
            self,
            settings: Settings,
            read_aligner: str) -> StrategyIndexAndAlignTumorNormal.main:

        indexer_class = self.INDEXER_DICT[read_aligner]
        aligner_class = self.ALIGNER_DICT[read_aligner]

        return StrategyIndexAndAlignTumorNormal(
            indexer=indexer_class(settings),
            aligner=aligner_class(settings)
        ).main
