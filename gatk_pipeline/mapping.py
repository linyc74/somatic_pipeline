from .constant import CMD_LINEBREAK
from .template import Processor, Settings


class BwaIndex(Processor):

    fna: str
    index: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fna: str) -> str:
        self.fna = fna

        self.set_index()
        self.execute()

        return self.index

    def set_index(self):
        self.index = f'{self.workdir}/bwa_index'

    def execute(self):
        stderr = f'{self.workdir}/bwa_index_stderr.log'

        cmd = CMD_LINEBREAK.join([
            'bwa index',
            f'-p {self.index}',
            self.fna,
            f'2> {stderr}',
        ])

        self.call(cmd)


class BwaMem(Processor):

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

        self.set_file_paths()
        self.mapping()
        self.sam_to_bam()
        self.sort_bam()

        return self.sorted_bam

    def set_file_paths(self):
        self.sam = f'{self.workdir}/{self.sample_name}.sam'
        self.bam = f'{self.workdir}/{self.sample_name}.bam'
        self.sorted_bam = f'{self.outdir}/{self.sample_name}_sorted.bam'

    def mapping(self):
        stderr = f'{self.workdir}/bwa_mem_{self.sample_name}_stderr.log'

        cmd = CMD_LINEBREAK.join([
            'bwa mem',
            f'-t {self.threads}',
            f'-o {self.sam}',
            self.index,
            self.fq1,
            self.fq2,
            f'2> {stderr}',
        ])

        self.call(cmd)

    def sam_to_bam(self):
        cmd = f'samtools view -b -h {self.sam} > {self.bam}'
        self.call(cmd)

    def sort_bam(self):
        cmd = f'samtools sort {self.bam} > {self.sorted_bam}'
        self.call(cmd)
