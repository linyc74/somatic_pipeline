import os
from typing import Tuple
from os.path import basename
from .template import Processor


class TrimGalore(Processor):

    QUALITY = 20
    LENGTH = 20
    MAX_N = 0
    CUTADAPT_TOTAL_CORES = 1
    # According to the help message of trim_galore, 2 cores for cutadapt -> actually up to 9 cores
    # Set it to 1 to avoid overloading the system when multiple jobs are running simultaneously

    fq1: str
    fq2: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    out_fq1: str
    out_fq2: str

    def main(
            self,
            fq1: str,
            fq2: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq1 = fq1
        self.fq2 = fq2
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.execute()
        self.move_fastqc_report()
        self.set_out_fq1()
        self.set_out_fq2()

        return self.out_fq1, self.out_fq2

    def execute(self):
        args = [
            'trim_galore',
            '--paired',
            f'--quality {self.QUALITY}',
            '--phred33',
            f'--cores {self.CUTADAPT_TOTAL_CORES}',
            f'--fastqc_args "--threads {self.threads}"',
            '--illumina',
            f'--length {self.LENGTH}',
            f'--max_n {self.MAX_N}',
            '--trim-n',
            '--gzip',
            f'--output_dir {self.workdir}'
        ]

        if self.clip_r1_5_prime > 0:
            args.append(f'--clip_R1 {self.clip_r1_5_prime}')

        if self.clip_r2_5_prime > 0:
            args.append(f'--clip_R2 {self.clip_r2_5_prime}')

        log = f'{self.outdir}/trim-galore.log'
        args += [
            self.fq1,
            self.fq2,
            f'1> {log} 2> {log}'
        ]

        self.call(self.CMD_LINEBREAK.join(args))

    def move_fastqc_report(self):
        dstdir = f'{self.outdir}/fastqc'
        os.makedirs(dstdir, exist_ok=True)
        for suffix in [
            'fastqc.html',
            'fastqc.zip',
            'trimming_report.txt'
        ]:
            self.call(f'mv {self.workdir}/*{suffix} {dstdir}/')

    def set_out_fq1(self):
        f = basename(self.fq1)
        f = self.__strip_file_extension(f)
        self.out_fq1 = f'{self.workdir}/{f}_val_1.fq.gz'

    def set_out_fq2(self):
        f = basename(self.fq2)
        f = self.__strip_file_extension(f)
        self.out_fq2 = f'{self.workdir}/{f}_val_2.fq.gz'

    def __strip_file_extension(self, f):
        for suffix in [
            '.fq',
            '.fq.gz',
            '.fastq',
            '.fastq.gz',
        ]:
            if f.endswith(suffix):
                f = f[:-len(suffix)]  # strip suffix
        return f
