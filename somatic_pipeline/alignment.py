from typing import Optional, Tuple
from abc import ABC, abstractmethod
from .template import Processor
from .constant import TUMOR, NORMAL


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
        dstdir = self.workdir if self.discard_bam else self.outdir
        self.sorted_bam = f'{dstdir}/{self.sample_name}-sorted.bam'
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
            f'-R "@RG\\tID:foo\\tSM:{self.sample_name}\\tPL:ILLUMINA"',  # needs literal '\t' for tab
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
            f'--rg-id foo',
            f'--rg "SM:{self.sample_name}"',
            f'--rg "PL:ILLUMINA"',
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
    discard_bam: bool

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
            normal_fq2: Optional[str],
            discard_bam: bool) -> Tuple[str, Optional[str]]:

        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.discard_bam = discard_bam

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
            sample_name=TUMOR,
            discard_bam=self.discard_bam)

    def align_normal(self):
        self.normal_bam = None if self.normal_fq1 is None \
            else self.aligner.main(
                index=self.index,
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                sample_name=NORMAL,
                discard_bam=self.discard_bam)


class FactoryIndexAndAlignTumorNormal(Processor):

    TO_INDEXER_CLASS = {
        'bwa': BwaIndexer,
        'bowtie2': Bowtie2Indexer,
    }
    TO_ALIGNER_CLASS = {
        'bwa': BwaAligner,
        'bowtie2': Bowtie2Aligner,
    }

    def get_strategy(
            self,
            read_aligner: str) -> StrategyIndexAndAlignTumorNormal:

        indexer_class = self.TO_INDEXER_CLASS[read_aligner]
        aligner_class = self.TO_ALIGNER_CLASS[read_aligner]

        return StrategyIndexAndAlignTumorNormal(
            indexer=indexer_class(self.settings),
            aligner=aligner_class(self.settings)
        )


class Alignment(Processor):

    read_aligner: str
    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    discard_bam: bool

    tumor_bam: str
    normal_bam: str

    def main(
            self,
            read_aligner: str,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str],
            discard_bam: bool) -> Tuple[str, str]:

        self.read_aligner = read_aligner
        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.discard_bam = discard_bam

        strategy = FactoryIndexAndAlignTumorNormal(self.settings).get_strategy(
            read_aligner=self.read_aligner)

        self.tumor_bam, self.normal_bam = strategy.main(
            ref_fa=self.ref_fa,
            tumor_fq1=self.tumor_fq1,
            tumor_fq2=self.tumor_fq2,
            normal_fq1=self.normal_fq1,
            normal_fq2=self.normal_fq2,
            discard_bam=self.discard_bam)

        return self.tumor_bam, self.normal_bam
