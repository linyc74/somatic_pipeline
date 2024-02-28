import gzip
from typing import Tuple, IO
from .tools import rev_comp
from .template import Processor


class TwistUMIConsensusMapping(Processor):
    """
    This class is an adaptation of the Twist UMI Consensus Mapping pipeline:
    https://www.twistbioscience.com/resources/guideguideline/processing-sequencing-data-utilizing-twist-unique-molecular-identifier-umi

    Commands in the PDF were buggy and had to be fixed.
    The fixed commands are in the `call_extremely_long_commands` method.

    This pipeline is not yet integrated into the main pipeline, due to its complexity.
    """

    ref_fa: str
    fq1: str
    fq2: str
    sample_name: str

    output_bam: str

    def main(
            self,
            ref_fa: str,
            fq1: str,
            fq2: str,
            sample_name: str) -> str:

        self.ref_fa = ref_fa
        self.fq1 = fq1
        self.fq2 = fq2
        self.sample_name = sample_name

        self.output_bam = f'{self.workdir}/{self.sample_name}_twist_umi_consensus_mapping.bam'
        self.call_extremely_long_commands()

        return self.output_bam

    def call_extremely_long_commands(self):
        cmds = [
            # (1) Convert Fastq to unaligned BAM
            [
                'picard FastqToSam',
                f'-O {self.workdir}/unaligned.bam',
                f'-F1 {self.fq1}',
                f'-F2 {self.fq2}',
                f'-SM {self.sample_name}',
                f'-LB {self.sample_name}',
                '-PU Unit1',
                '-PL Illumina',
            ],
            # (2) Extract UMI bases into unaligned BAM tag
            [
                'fgbio ExtractUmisFromBam',
                f'--input={self.workdir}/unaligned.bam',
                f'--output={self.workdir}/unaligned_umi_extracted.bam',
                '--read-structure=5M2S+T 5M2S+T',
                '--molecular-index-tags=ZA ZB',
                '--single-tag=RX',
            ],
            # (3) Convert unaligned BAM to Fastq
            [
                'picard SamToFastq',
                f'--INPUT {self.workdir}/unaligned_umi_extracted.bam',
                f'--FASTQ {self.workdir}/umi_extracted.fastq',
                '--INTERLEAVE true',
            ],
            # (4) Align Fastq
            [
                f'bwa mem -p -t {self.threads}',
                f'{self.ref_fa}',
                f'{self.workdir}/umi_extracted.fastq',
                '|',
                f'samtools sort -@ {self.threads}',
                f'-o {self.workdir}/aligned_umi_extracted.bam',
            ],
            # (5) Merge aligned BAM with unaligned BAM containing UMI tags
            [
                'picard MergeBamAlignment',
                f'--UNMAPPED_BAM {self.workdir}/unaligned_umi_extracted.bam',
                f'--ALIGNED_BAM {self.workdir}/aligned_umi_extracted.bam',
                f'--OUTPUT {self.workdir}/aligned_tag_umi.bam',
                f'--REFERENCE_SEQUENCE {self.ref_fa}',
                '--CLIP_ADAPTERS false',
                '--VALIDATION_STRINGENCY SILENT',
                '--CREATE_INDEX true',
                '--EXPECTED_ORIENTATIONS FR',
                '--MAX_GAPS -1',
                '--SORT_ORDER coordinate',
                '--ALIGNER_PROPER_PAIR_FLAGS false',
            ],
            # (6) Optional: Collect Hs metrics prior to consensus calling
            # (7) Group reads by UMI
            [
                'fgbio GroupReadsByUmi',
                '--strategy=Paired',
                f'--input={self.workdir}/aligned_tag_umi.bam',
                f'--output={self.workdir}/grouped_by_umi.bam',
                '--raw-tag=RX',
                '--min-map-q=10',
                '--edits=1',
            ],
            # (8) Call consensus reads
            [
                'fgbio CallDuplexConsensusReads',
                f'--input={self.workdir}/grouped_by_umi.bam',
                f'--output={self.workdir}/unaligned_consensus.bam',
                '--error-rate-pre-umi=45',
                '--error-rate-post-umi=30',
                '--min-input-base-quality=30',
                '--min-reads 2 1 1',
            ],
            # (9) Align duplex consensus reads
            [
                'picard SamToFastq',
                f'--INPUT {self.workdir}/unaligned_consensus.bam',
                f'--FASTQ {self.workdir}/consensus.fastq',
                '--INTERLEAVE true',
            ],
            [
                f'bwa mem -p -t {self.threads}',
                f'{self.ref_fa}',
                f'-o {self.workdir}/aligned_consensus.bam',
                f'{self.workdir}/consensus.fastq',
            ],

            # The following SortSam step needs to be added to prevent MergeBamAlignment error:
            # ```Aligned record iterator (...) is behind the unmapped reads (...)```

            # As stated in https://sourceforge.net/p/samtools/mailman/samtools-help/thread/d720f26439c8e4af960b0d3c37ba6918.squirrel@webmail.erasmusmc.nl/:
            # ```
            #    Note that there is a long-standing discrepancy between the way samtools and
            #    Picard define queryname order, so it could be that you need to sort your
            #    unmapped BAM with Picard SortSam.  If the header says that it is queryname
            #    sorted, you can run ValidateSam on it to confirm that it is queryname
            #    sorted.
            # ```
            [
                'picard SortSam',
                f'--INPUT {self.workdir}/unaligned_consensus.bam',
                f'--OUTPUT {self.workdir}/unaligned_consensus_queryname_sorted.bam',
                '-SORT_ORDER queryname',
            ],
            # (10) Merge unaligned consensus BAM and aligned consensus BAM to retain UMI tag metadata
            [
                'picard MergeBamAlignment',
                f'--UNMAPPED_BAM {self.workdir}/unaligned_consensus_queryname_sorted.bam',
                f'--ALIGNED_BAM {self.workdir}/aligned_consensus.bam',
                f'--OUTPUT {self.workdir}/merged_no_read_group_consensus.bam',
                f'--REFERENCE_SEQUENCE {self.ref_fa}',
                '--CLIP_ADAPTERS false',
                '--VALIDATION_STRINGENCY SILENT',
                '--CREATE_INDEX true',
                '--EXPECTED_ORIENTATIONS FR',
                '--MAX_GAPS -1',
                '--SORT_ORDER coordinate',
                '--ALIGNER_PROPER_PAIR_FLAGS false',
                '--ATTRIBUTES_TO_RETAIN X0',
                '--ATTRIBUTES_TO_RETAIN ZS',
                '--ATTRIBUTES_TO_RETAIN ZI',
                '--ATTRIBUTES_TO_RETAIN ZM',
                '--ATTRIBUTES_TO_RETAIN ZC',
                '--ATTRIBUTES_TO_RETAIN ZN',
                '--ATTRIBUTES_TO_RETAIN ad',
                '--ATTRIBUTES_TO_RETAIN bd',
                '--ATTRIBUTES_TO_RETAIN cd',
                '--ATTRIBUTES_TO_RETAIN ae',
                '--ATTRIBUTES_TO_RETAIN be',
                '--ATTRIBUTES_TO_RETAIN ce',
            ],
            [
                'picard AddOrReplaceReadGroups',
                f'--INPUT {self.workdir}/merged_no_read_group_consensus.bam',
                f'--OUTPUT {self.output_bam}',
                f'--RGID {self.sample_name}',
                f'--RGLB {self.sample_name}',
                '--RGPL Illumina',
                f'--RGSM {self.sample_name}',
                '--RGPU NA',
            ],
        ]
        for cmd in cmds:
            log = f'{self.outdir}/twist_mapping.log'
            cmd.append(f'1>> {log} 2>> {log}')
            self.call(self.CMD_LINEBREAK.join(cmd))


class RemoveTwistUMI(Processor):

    UMI_LENGTH = 7

    fq1: str
    fq2: str
    gz: bool

    out_fq1: str
    out_fq2: str

    fq1_reader: IO
    fq2_reader: IO
    fq1_writer: IO
    fq2_writer: IO

    def main(self, fq1: str, fq2: str, gz: bool) -> Tuple[str, str]:
        self.fq1 = fq1
        self.fq2 = fq2
        self.gz = gz

        self.open_files()

        self.logger.info(f'Removing Twist UMI and universal adapter from:\n{fq1} -> {self.out_fq1}\n{fq2} -> {self.out_fq2}')

        while True:
            header1 = self.fq1_reader.readline().strip()
            header2 = self.fq2_reader.readline().strip()

            if header1 == '':  # end of file
                break

            assert header1.split()[0] == header2.split()[0]

            seq1 = self.fq1_reader.readline().strip()[self.UMI_LENGTH:]
            seq2 = self.fq2_reader.readline().strip()[self.UMI_LENGTH:]
            new_seq2 = strip_mate_3prime_umi(read=seq1, mate=seq2)
            new_seq1 = strip_mate_3prime_umi(read=seq2, mate=seq1)

            self.fq1_reader.readline()  # skip the '+' line
            self.fq2_reader.readline()

            qual1 = self.fq1_reader.readline().strip()
            qual2 = self.fq2_reader.readline().strip()
            u = self.UMI_LENGTH
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
        self.gzip_output()

        return self.out_fq1, self.out_fq2

    def open_files(self):
        self.out_fq1 = f'{self.workdir}/twist_umi_removed_R1.fastq'
        self.out_fq2 = f'{self.workdir}/twist_umi_removed_R2.fastq'
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
            self.call(self.out_fq1)
            self.call(self.out_fq2)
            self.out_fq1 += '.gz'
            self.out_fq2 += '.gz'


def strip_mate_3prime_umi(read: str, mate: str) -> str:
    mate_rc = rev_comp(mate)
    pos = mate_rc.find(read[:15])  # 4^15 = 1,073,741,824 should be specific enough
    return mate[:-pos] if pos > 0 else mate
