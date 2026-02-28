# VEP Installation

Stand-alone installation of VEP version 115.2

```bash
# use the perl in the conda environment
conda install -c conda-forge perl-list-moreutils
conda install -c conda-forge perl-libwww-perl
conda install -c bioconda perl-dbi
conda install -c bioconda perl-bio-db-hts

cd ~/opt
wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/115.2.tar.gz
tar xzf 115.2.tar.gz
rm 115.2.tar.gz
cd ensembl-vep-release-115.2

# include "--NO_UPDATE" so that the installation process will not be disrupted by update check
# disruption of installation will result in missing perl modules (e.g. Bio/EnsEMBL/Registry.pm) and plugins
perl INSTALL.pl --AUTO ap --PLUGINS all --NO_HTSLIB --NO_UPDATE

export PATH=$PATH:$HOME/opt/ensembl-vep-release-115.2  # in .bashrc
```
