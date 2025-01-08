import os
from .template import Processor


class PCGR(Processor):

    GENOME_ASSEMBLY = 'grch38'
    SEQUENCING_ASSAY = 'WES'

    vcf: str
    pcgr_ref_data_tgz: str
    pcgr_vep_tar_gz: str
    vep_buffer_size: int
    pcgr_tumor_site: int
    pcgr_tmb_target_size_mb: int
    pcgr_tmb_display: str

    refdata_dir: str
    vep_dir: str

    def main(
            self,
            vcf: str,
            pcgr_ref_data_tgz: str,
            pcgr_vep_tar_gz: str,
            vep_buffer_size: int,
            pcgr_tumor_site: int,
            pcgr_tmb_target_size_mb: int,
            pcgr_tmb_display: str):

        self.vcf = vcf
        self.pcgr_ref_data_tgz = pcgr_ref_data_tgz
        self.pcgr_vep_tar_gz = pcgr_vep_tar_gz
        self.vep_buffer_size = vep_buffer_size
        self.pcgr_tumor_site = pcgr_tumor_site
        self.pcgr_tmb_target_size_mb = pcgr_tmb_target_size_mb
        self.pcgr_tmb_display = pcgr_tmb_display

        self.untar_reference_data()
        self.index_vcf()
        self.run_pcgr()

    def untar_reference_data(self):
        self.call(f'tar xzf {self.pcgr_ref_data_tgz} -C {self.workdir}')  # creates data/ in workdir
        self.refdata_dir = self.workdir
        # workdir/data/grch38/...
        # ^^^^^^^  the PCGR refdata dir should be at this level, which contains the 'data/' dir

        if os.path.isfile(self.pcgr_vep_tar_gz):
            self.call(f'tar xzf {self.pcgr_vep_tar_gz} -C {self.workdir}')  # creates homo_sapiens/ in workdir
            self.vep_dir = self.workdir
            # workdir/homo_sapiens/112_GRCh38/...
            # ^^^^^^^  the VEP cache dir should be at this level, which contains the 'homo_sapiens/' dir

        else:  # already extracted, with `homo_sapiens` placed in the `vep_cache` dir
            self.vep_dir = self.pcgr_vep_tar_gz

    def index_vcf(self):
        if not os.path.exists(f'{self.vcf}.tbi'):
            self.call(f'tabix --preset vcf {self.vcf}')

    def run_pcgr(self):
        log = f'{self.outdir}/pcgr.log'
        cmd = self.CMD_LINEBREAK.join([
            'pcgr',
            f'--input_vcf {self.vcf}',
            f'--vep_dir {self.vep_dir}',
            f'--refdata_dir {self.refdata_dir}',
            f'--output_dir {self.outdir}/pcgr',
            f'--genome_assembly {self.GENOME_ASSEMBLY}',
            f'--assay {self.SEQUENCING_ASSAY}',
            f'--sample_id {os.path.basename(self.outdir)}',

            f'--vcfanno_n_proc {self.threads}',
            f'--vep_n_forks {self.threads}',
            f'--vep_buffer_size {self.vep_buffer_size}',

            f'--effective_target_size_mb {self.pcgr_tmb_target_size_mb}',
            f'--estimate_tmb',
            f'--tmb_display {self.pcgr_tmb_display}',

            f'--estimate_signatures',  # Estimate relative contributions of reference mutational signatures in query sample (re-fitting), default: False

            f'--tumor_site {self.pcgr_tumor_site}',

            f'1>> {log} 2>> {log}',
        ])
        self.call(cmd)
