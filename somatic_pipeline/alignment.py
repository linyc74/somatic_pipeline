from typing import Optional, Tuple
from .constant import TUMOR, NORMAL
from .template import Processor, Settings


class Indexer(Processor):

    fna: str
    index: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fna: str) -> str:
        self.fna = fna
        self.execute()
        return self.index

    def execute(self):
        pass


class BwaIndexer(Indexer):

    def execute(self):
        self.index = f'{self.workdir}/bwa-index'
        log = f'{self.outdir}/bwa-index.log'
        cmd = self.CMD_LINEBREAK.join([
            'bwa index',
            f'-p {self.index}',
            self.fna,
            f'2> {log}',
        ])
        self.call(cmd)


class Bowtei2Indexer(Indexer):

    def execute(self):
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


class Aligner(Processor):

    index: str
    fq1: str
    fq2: str
    sample_name: str

    sam: str
    bam: str
    sorted_bam: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

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

    def align(self):
        pass

    def sam_to_bam(self):
        self.bam = f'{self.workdir}/{self.sample_name}.bam'
        cmd = f'samtools view -b -h {self.sam} > {self.bam}'
        self.call(cmd)

    def sort_bam(self):
        self.sorted_bam = f'{self.outdir}/{self.sample_name}-sorted.bam'
        cmd = f'samtools sort {self.bam} > {self.sorted_bam}'
        self.call(cmd)


class BwaAligner(Aligner):

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
            f'--threads {self.threads}',
            f'2> {log}',
        ])
        self.call(cmd)


class Alignment(Processor):

    INDEXER_DICT = {
        'bwa': BwaIndexer,
        'bowtie2': Bowtei2Indexer,
    }
    ALIGNER_DICT = {
        'bwa': BwaAligner,
        'bowtie2': Bowtie2Aligner,
    }

    read_aligner: str
    ref_fa: str
    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]

    index: str
    aligner: Aligner
    tumor_bam: str
    normal_bam: Optional[str]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            read_aligner: str,
            ref_fa: str,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str]) -> Tuple[str, Optional[str]]:

        self.read_aligner = read_aligner
        self.ref_fa = ref_fa
        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2

        self.index_genome()
        self.set_aligner()
        self.align_tumor()
        self.align_normal()

        return self.tumor_bam, self.normal_bam

    def index_genome(self):
        indexer = self.INDEXER_DICT[self.read_aligner](self.settings)
        self.index = indexer.main(
            fna=self.ref_fa)

    def set_aligner(self):
        self.aligner = self.ALIGNER_DICT[self.read_aligner](self.settings)

    def align_tumor(self):
        self.tumor_bam = self.aligner.main(
            index=self.index,
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2,
            sample_name=TUMOR)

    def align_normal(self):
        if self.normal_fq1 is None:
            self.normal_bam = None
        else:
            self.normal_bam = self.aligner.main(
                index=self.index,
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                sample_name=NORMAL)
