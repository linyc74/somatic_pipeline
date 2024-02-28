import os
import gzip
from typing import Tuple, Optional, IO, List
from os.path import basename
from .template import Processor
from .tools import rev_comp, edit_fpath


class Trimming(Processor):

    tumor_fq1: str
    tumor_fq2: str
    normal_fq1: Optional[str]
    normal_fq2: Optional[str]
    umi_length: int
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    files_to_be_removed: List[str]

    def main(
            self,
            tumor_fq1: str,
            tumor_fq2: str,
            normal_fq1: Optional[str],
            normal_fq2: Optional[str],
            umi_length: int,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str, Optional[str], Optional[str]]:

        self.tumor_fq1 = tumor_fq1
        self.tumor_fq2 = tumor_fq2
        self.normal_fq1 = normal_fq1
        self.normal_fq2 = normal_fq2
        self.umi_length = umi_length
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.files_to_be_removed = []

        if self.umi_length > 0:
            self.tumor_fq1, self.tumor_fq2 = RemoveUmiAndAdapter(self.settings).main(
                fq1=self.tumor_fq1,
                fq2=self.tumor_fq2,
                umi_length=self.umi_length,
                gz=False)
            self.files_to_be_removed += [self.tumor_fq1, self.tumor_fq2]

            if self.normal_fq1 is not None:
                self.normal_fq1, self.normal_fq2 = RemoveUmiAndAdapter(self.settings).main(
                    fq1=self.normal_fq1,
                    fq2=self.normal_fq2,
                    umi_length=self.umi_length,
                    gz=False)
                self.files_to_be_removed += [self.normal_fq1, self.normal_fq2]

        self.tumor_fq1, self.tumor_fq2 = TrimGalore(self.settings).main(
            fq1=self.tumor_fq1,
            fq2=self.tumor_fq2,
            clip_r1_5_prime=self.clip_r1_5_prime,
            clip_r2_5_prime=self.clip_r2_5_prime)

        if self.normal_fq1 is not None:
            self.normal_fq1, self.normal_fq2 = TrimGalore(self.settings).main(
                fq1=self.normal_fq1,
                fq2=self.normal_fq2,
                clip_r1_5_prime=self.clip_r1_5_prime,
                clip_r2_5_prime=self.clip_r2_5_prime)

        for file in self.files_to_be_removed:
            self.call(f'rm {file}')

        return self.tumor_fq1, self.tumor_fq2, self.normal_fq1, self.normal_fq2


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
        f = self.__strip_fastq_extension(f)
        self.out_fq1 = f'{self.workdir}/{f}_val_1.fq.gz'

    def set_out_fq2(self):
        f = basename(self.fq2)
        f = self.__strip_fastq_extension(f)
        self.out_fq2 = f'{self.workdir}/{f}_val_2.fq.gz'

    def __strip_fastq_extension(self, f):
        for suffix in [
            '.fq',
            '.fq.gz',
            '.fastq',
            '.fastq.gz',
        ]:
            if f.endswith(suffix):
                f = f[:-len(suffix)]  # strip suffix
        return f


class RemoveUmiAndAdapter(Processor):

    fq1: str
    fq2: str
    umi_length: int
    gz: bool

    out_fq1: str
    out_fq2: str

    fq1_reader: IO
    fq2_reader: IO
    fq1_writer: IO
    fq2_writer: IO

    def main(self, fq1: str, fq2: str, umi_length: int, gz: bool) -> Tuple[str, str]:

        self.fq1 = fq1
        self.fq2 = fq2
        self.umi_length = umi_length
        self.gz = gz

        self.open_files()

        self.logger.info(f'Start removing UMI ({self.umi_length} bp) and universal adapters of "{self.fq1}" and "{self.fq2}"...')

        total_1, total_2, remain_1, remain_2 = 0, 0, 0, 0
        while True:
            header1 = self.fq1_reader.readline().strip()
            header2 = self.fq2_reader.readline().strip()

            if header1 == '':  # end of file
                break

            assert header1.split()[0] == header2.split()[0]

            seq1 = self.fq1_reader.readline().strip()
            seq2 = self.fq2_reader.readline().strip()
            total_1 += len(seq1)
            total_2 += len(seq2)

            seq1 = seq1[self.umi_length:]  # 5' clip
            seq2 = seq2[self.umi_length:]
            new_seq2 = strip_mate_3prime_umi(read=seq1, mate=seq2)
            new_seq1 = strip_mate_3prime_umi(read=seq2, mate=seq1)
            remain_1 += len(new_seq1)
            remain_2 += len(new_seq2)

            self.fq1_reader.readline()  # skip the '+' line
            self.fq2_reader.readline()

            qual1 = self.fq1_reader.readline().strip()
            qual2 = self.fq2_reader.readline().strip()
            u = self.umi_length
            qual1 = qual1[u:u+len(new_seq1)]
            qual2 = qual2[u:u+len(new_seq2)]

            assert len(new_seq1) == len(qual1)
            assert len(new_seq2) == len(qual2)

            self.fq1_writer.write(header1 + '\n')
            self.fq2_writer.write(header2 + '\n')
            self.fq1_writer.write(new_seq1 + '\n')
            self.fq2_writer.write(new_seq2 + '\n')
            self.fq1_writer.write('+\n')
            self.fq2_writer.write('+\n')
            self.fq1_writer.write(qual1 + '\n')
            self.fq2_writer.write(qual2 + '\n')

        self.close_files()

        self.logger.info(f'''\
{self.fq1} ({total_1:,} bp) -> ({remain_1:,} bp = {remain_1/total_1*100:.2f}%) {self.out_fq1}
{self.fq2} ({total_2:,} bp) -> ({remain_2:,} bp = {remain_2/total_2*100:.2f}%) {self.out_fq2}''')

        self.gzip_output()

        return self.out_fq1, self.out_fq2

    def open_files(self):

        self.out_fq1 = edit_fpath(
            fpath=self.fq1,
            old_suffix=get_fastq_ext(self.fq1),
            new_suffix='_umi_adapter_removed.fastq',
            dstdir=self.workdir)

        self.out_fq2 = edit_fpath(
            fpath=self.fq2,
            old_suffix=get_fastq_ext(self.fq2),
            new_suffix='_umi_adapter_removed.fastq',
            dstdir=self.workdir)

        self.fq1_reader = gzip.open(self.fq1, 'rt') if self.fq1.endswith('.gz') else open(self.fq1, 'r')
        self.fq2_reader = gzip.open(self.fq2, 'rt') if self.fq2.endswith('.gz') else open(self.fq2, 'r')
        self.fq1_writer = open(self.out_fq1, 'w')
        self.fq2_writer = open(self.out_fq2, 'w')

    def close_files(self):
        self.fq1_reader.close()
        self.fq2_reader.close()
        self.fq1_writer.close()
        self.fq2_writer.close()

    def gzip_output(self):
        if self.gz:
            self.call(f'gzip {self.out_fq1}')
            self.call(f'gzip {self.out_fq2}')
            self.out_fq1 += '.gz'
            self.out_fq2 += '.gz'


def strip_mate_3prime_umi(read: str, mate: str) -> str:
    mate_rc = rev_comp(mate)
    pos = mate_rc.find(read[:15])  # 4^15 = 1,073,741,824 should be specific enough
    return mate[:-pos] if pos > 0 else mate


def get_fastq_ext(fq: str) -> str:
    for suffix in [
        '.fq',
        '.fq.gz',
        '.fastq',
        '.fastq.gz',
    ]:
        if fq.endswith(suffix):
            return suffix
    return ''
