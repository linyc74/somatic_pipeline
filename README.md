# Somatic Variant Calling Pipeline

## Usage

```bash
git clone https://github.com/linyc74/somatic_pipeline.git
```

Tumor-normal paired mode

```bash
python somatic_pipeline \
  -r reference_genome.fa \
  -1 tumor.1.fq.gz \
  -2 tumor.2.fq.gz \
  -3 normal.1.fq.gz \
  -4 normal.1.fq.gz
```

Tumor-only mode

```bash
python somatic_pipeline \
  -r reference_genome.fa \
  -1 tumor.1.fq.gz \
  -2 tumor.2.fq.gz
```

## Environment

Dependencies:
- `trim-galore`
- `bwa`
- `bowtie2`
- `samtools`
- `bcftools`
- `bedtools`
- `GATK`
- `MuSE`
- `VarScan`
- `LoFreq`
- `SomaticSniper`
- `VarDictJava`
- `SnpEff`
- `VEP`
- `CNVkit`
- `vcf2maf`
- `pandas`
