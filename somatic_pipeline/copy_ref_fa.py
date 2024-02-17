from os.path import basename, exists
from .template import Processor


class CopyRefFa(Processor):

    ref_fa: str
    copied_ref_fa: str

    def main(self, ref_fa: str) -> str:
        self.ref_fa = ref_fa

        self.unzip_or_copy()
        self.copy_bwa_index_files()

        return self.copied_ref_fa

    def unzip_or_copy(self):
        fname = basename(self.ref_fa).rstrip('.gz')
        dst = f'{self.workdir}/{fname}'

        if self.ref_fa.endswith('.gz'):
            self.call(f'gzip --decompress --stdout {self.ref_fa} > {dst}')
        else:
            self.call(f'cp {self.ref_fa} {dst}')

        self.copied_ref_fa = dst

    def copy_bwa_index_files(self):
        fname = self.ref_fa.rstrip('.gz')
        for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
            if exists(f'{fname}{ext}'):
                cmd = f'cp {fname}{ext} {self.workdir}/bwa-index{ext}'
                self.call(cmd)
