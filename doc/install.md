# Installation

Use PCGR as the base environment for everything.

```bash
export PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v2.2.5/conda/env/lock/"
export PLATFORM="linux"
conda create --name pcgr --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock
conda create --name pcgrr --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock
```

PCGR by default includes the following: `samtools`, `bcftools`, `vcf2maf`, `bedtools`.
Install the remaining. These versions have been tested thoroughly.
In particular, `muse=1.0` and `snpsift=5.2` are not the latest versions as their later versions cause bugs or dependency issues.

```bash
conda install -c bioconda \
  bwa=0.7.18 \
  bowtie2=2.5.4 \
  gatk4=4.3.0 \
  lofreq=2.1.5 \
  varscan=2.4.6 \
  somatic-sniper \
  muse=1.0 \
  snpsift=5.2
```

Re-install latest version of `bcftools` to avoid segmentation fault.

```bash
# only remove bcftools v1.21 without also removing lofreq
conda remove --force-remove bcftools

# https://www.htslib.org/download/
cd ~/opt
wget https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2
bzip2 -d bcftools-1.23.tar.bz2
tar xf bcftools-1.23.tar
rm bcftools-1.23.tar
cd bcftools-1.23
./configure
make

export PATH=$PATH:$HOME/opt/bcftools-1.23  # in .bashrc
```

Trim-Galore cannot be installed directly in the PCGR environment. Directly download from source.

```bash
pip install cutadapt==5.2
conda install -c bioconda fastqc=0.12.1

cd ~/opt
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
rm trim_galore.tar.gz

export PATH=$PATH:$HOME/opt/TrimGalore-0.6.10  # in .bashrc
```

Install VarDict from source.

```bash
sudo apt install dos2unix

cd ~/opt
git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git
cd VarDictJava
dos2unix gradlew
./gradlew distTar
tar -xf build/distributions/VarDict-1.8.3.tar
cp -r VarDict-1.8.3 ~/opt/
cd ..
rm -rf VarDictJava

dos2unix ~/opt/VarDict-1.8.3/bin/*.pl
dos2unix ~/opt/VarDict-1.8.3/bin/*.R

export PATH=$PATH:$HOME/opt/VarDict-1.8.3/bin  # in .bashrc
```

Install MSIsensor from source to avoid segmentation fault.

```bash
cd ~/opt
wget https://github.com/ding-lab/msisensor/releases/download/0.6/msisensor.linux
chmod 755 msisensor.linux
mv msisensor.linux msisensor

export PATH=$PATH:$HOME/opt  # in .bashrc
```

Install CNVkit.

```bash
pip install cnvkit==0.9.13

conda install -n pcgrr -c bioconda bioconductor-dnacopy

# link the R/Rscript in pcgrr to pcgr environment
ENVS=$HOME/anaconda3/envs
ln -s $ENVS/pcgrr/bin/R $ENVS/pcgr/bin/R
ln -s $ENVS/pcgrr/bin/Rscript $ENVS/pcgr/bin/Rscript
```
