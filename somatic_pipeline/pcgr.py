import os
from .template import Processor


class PCGR(Processor):

    GENOME_ASSEMBLY = 'grch38'
    SEQUENCING_ASSAY = 'WES'
    VEP_BUFFER_SIZE = 100
    EFFECTIVE_TARGET_SIZE_MB = 34  # Effective target size in Mb (potentially limited by read depth) of sequencing assay (for TMB analysis) (default: 34 (WES/WGS))
    TMB_DISPLAY = 'coding_and_silent'  # type of TMB measure to show in report: coding_and_silent, coding_non_silent, missense_only
    TUMOR_SITE = 12  # head and neck

    vcf: str
    pcgr_ref_data_tgz: str
    pcgr_vep_tar_gz: str

    refdata_dir: str
    vep_dir: str

    def main(
            self,
            vcf: str,
            pcgr_ref_data_tgz: str,
            pcgr_vep_tar_gz: str):

        self.vcf = vcf
        self.pcgr_ref_data_tgz = pcgr_ref_data_tgz
        self.pcgr_vep_tar_gz = pcgr_vep_tar_gz

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
            f'--vep_buffer_size {self.VEP_BUFFER_SIZE}',

            f'--effective_target_size_mb {self.EFFECTIVE_TARGET_SIZE_MB}',
            f'--estimate_tmb',
            f'--tmb_display {self.TMB_DISPLAY}',

            f'--estimate_signatures',  # Estimate relative contributions of reference mutational signatures in query sample (re-fitting), default: False

            f'--tumor_site {self.TUMOR_SITE}',
        ])
        self.call(cmd)
