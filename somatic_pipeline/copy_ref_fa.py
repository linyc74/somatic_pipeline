from os.path import basename
from .template import Processor


class CopyRefFa(Processor):

    ref_fa: str
    copied_ref_fa: str

    def main(self, ref_fa: str) -> str:
        self.ref_fa = ref_fa

        self.set_copied_fa()
        self.unzip_or_copy()

        return self.copied_ref_fa

    def set_copied_fa(self):
        fname = basename(self.ref_fa)
        if self.ref_fa.endswith('.gz'):
            fname = fname[:-3]
        self.copied_ref_fa = f'{self.workdir}/{fname}'

    def unzip_or_copy(self):
        if self.ref_fa.endswith('.gz'):
            cmd = f'gzip --decompress --stdout {self.ref_fa} > {self.copied_ref_fa}'
        else:
            cmd = f'cp {self.ref_fa} {self.copied_ref_fa}'
        self.call(cmd)
