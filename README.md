# GATK Pipeline

**Custom-built GATK pipeline**

## Usage

```bash
git clone https://github.com/linyc74/gatk_pipeline.git
```

Tumor-normal paired mode

```bash
python gatk_pipeline \
  -r reference_genome.fa \
  -1 tumor.1.fq.gz \
  -2 tumor.2.fq.gz \
  -3 normal.1.fq.gz \
  -4 normal.1.fq.gz
```

Tumor-only mode

```bash
python gatk_pipeline \
  -r reference_genome.fa \
  -1 tumor.1.fq.gz \
  -2 tumor.2.fq.gz
```

## Environment

Create a conda environment and install the following packages:

```bash
conda create --name gatk
conda activate gatk
conda install -c bioconda trim-galore bwa samtools gatk4
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
