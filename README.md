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

Required packages:
- `trim-galore`
- `bwa`
- `bowtie2`
- `samtools`
- `bcftools`
- `GATK`
- `MuSE`
- `VarScan`
- `SnpEff`
- `pandas`

Create a conda environment `somatic` and install packages:

```bash
conda create --name somatic
conda activate somatic
conda install -c bioconda trim-galore bwa bowtie2 samtools bcftools gatk4 muse varscan
pip install pandas
```

Install SnpEff:

```bash
mkdir ~/opt
cd ~/opt

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

# Add snpEff to PATH variable
export PATH=$PATH:$HOME/opt/snpEff/exec
```
