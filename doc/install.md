## Conda

```bash
conda create -n somatic python=3.10
conda activate somatic

conda install -c anaconda pandas=2.1.4
conda install -c conda-forge tbb=2020.2 unzip=6.0
conda install -c bioconda \
trim-galore=0.6.6 \
bwa=0.7.18 \
samtools=1.21 \
gatk4=4.6.1.0 \
bowtie2=2.5.4 \
muse=1.0 \
varscan=2.3.7 \
bcftools=1.11 \
vcf2maf=1.6.21 \
bedtools=2.30.0 \
somatic-sniper=1.0.5.0

conda clean --all --yes
```

## LoFreq

```bash
cd ~/opt
wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz
tar -xzf lofreq_star-2.1.5_linux-x86-64.tgz
rm lofreq_star-2.1.5_linux-x86-64.tgz
```

In `.bashrc` add:
```bash
export PATH=$PATH:$HOME/opt/lofreq_star-2.1.5_linux-x86-64/bin
```

## VarDict

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
```

In `.bashrc` add:
```bash
export PATH=$PATH:$HOME/opt/VarDict-1.8.3/bin
```

## VEP

```bash
RUN apt-get install -y \
    build-essential \
    perl \
    cpanminus \
    curl \
    libwww-perl \
 && cpanm DBI \
 && cpanm Try::Tiny

cd ~/opt
wget https://github.com/Ensembl/ensembl-vep/archive/release/106.zip
unzip 106.zip
rm 106.zip
cd ensembl-vep-release-106

# include "--NO_UPDATE" so that the installation process will not be disrupted by update check
# disruption of installation will result in missing perl modules (e.g. Bio/EnsEMBL/Registry.pm) and plugins
perl INSTALL.pl --AUTO ap --PLUGINS all --NO_HTSLIB --NO_UPDATE
```

In `.bashrc` add:
```bash
export PATH=$PATH:$HOME/opt/ensembl-vep-release-106
```

## PCGR

```bash
PCGR_VERSION="2.1.2"
PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock/"
PLATFORM="linux"
conda create --name pcgr --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock
conda create --name pcgrr --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock
```

In `.bashrc` add the following line to make `pcgr` executable:
```bash
export PATH=$PATH:$HOME/anaconda3/envs/pcgr/bin
```

Make `pcgrr.R` available in the current `somatic` environment:
```bash
ENVS=$HOME/anaconda3/envs
cp $ENVS/pcgr/bin/pcgrr.R $ENVS/somatic/bin/
```
