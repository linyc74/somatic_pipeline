from os.path import exists
from .tools import edit_fpath
from .template import Processor


class SamtoolsIndexFa(Processor):

    fa: str

    def main(self, fa: str):
        self.fa = fa
        if self.index_file_already_exist():
            return
        self.execute()

    def index_file_already_exist(self) -> bool:
        return exists(f'{self.fa}.fai')

    def execute(self):
        self.call(f'samtools faidx {self.fa}')


class SamtoolsIndexBam(Processor):

    bam: str

    def main(self, bam: str):
        self.bam = bam
        if self.index_file_already_exist():
            return
        self.execute()

    def index_file_already_exist(self) -> bool:
        return exists(f'{self.bam}.bai')

    def execute(self):
        self.call(f'samtools index {self.bam}')


class GATKCreateSequenceDictionary(Processor):

    ref_fa: str

    def main(self, ref_fa: str):
        self.ref_fa = ref_fa
        if self.dict_file_already_exist():
            return
        self.execute()

    def dict_file_already_exist(self) -> bool:
        dict_file = edit_fpath(
            fpath=self.ref_fa,
            old_suffix='.fa',
            new_suffix='.dict')
        return exists(dict_file)

    def execute(self):
        log = f'{self.outdir}/gatk-CreateSequenceDictionary.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk CreateSequenceDictionary',
            f'-R {self.ref_fa}',
            f'1>> {log}',
            f'2>> {log}',
        ])
        self.call(cmd)


class GATKIndexVcf(Processor):

    vcf: str

    def main(self, vcf: str):
        self.vcf = vcf
        if self.index_file_already_exist():
            self.logger.debug(f'Index file for {self.vcf} exists, skip indexing')
            return
        self.execute()

    def index_file_already_exist(self) -> bool:
        return exists(f'{self.vcf}.idx') or exists(f'{self.vcf}.tbi')

    def execute(self):
        log = f'{self.outdir}/gatk-IndexFeatureFile.log'
        cmd = self.CMD_LINEBREAK.join([
            'gatk IndexFeatureFile',
            f'--input {self.vcf}',
            f'1>> {log}',
            f'2>> {log}',
        ])
        self.call(cmd)


class BgzipIndex(Processor):

    vcf: str
    keep: bool

    vcf_gz: str

    def main(self, vcf: str, keep: bool = True) -> str:
        self.vcf = vcf
        self.keep = keep

        self.block_gzip()
        self.tabix_index()

        return self.vcf_gz

    def block_gzip(self):
        self.vcf_gz = self.vcf + '.gz'
        self.call(f'bgzip --threads {self.threads} --stdout {self.vcf} > {self.vcf_gz}')
        if not self.keep:
            self.call(f'rm {self.vcf}')

    def tabix_index(self):
        self.call(f'tabix --preset vcf {self.vcf_gz}')
