from .template import Processor


class TwistUMIConsensusMapping(Processor):
    """
    This class is an adaptation of the Twist UMI Consensus Mapping pipeline:
    https://www.twistbioscience.com/resources/guideguideline/processing-sequencing-data-utilizing-twist-unique-molecular-identifier-umi

    Commands in the PDF were buggy and had to be fixed.
    The fixed commands are in the `call_extremely_long_commands` method.
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
