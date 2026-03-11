# Perl & VEP (Painful) Setup

Ensembl Variant Effect Predictor (VEP) was built on Perl.
Installing VEP requires a bunch of Perl dependencies,
which is much harder to install compared to Python packages.

## Basics

Perl packages, unlike Conda, are not pre-built, so they need to be compiled.

First, install the "correct" Perl version from conda:

```bash
conda remove perl
conda install -c conda-forge perl=5.26.2=h470a237_0
```

This specific build `h470a237_0`, does not contain hard-coded `gcc` path in its config file.
Hard-coding path is the developer's bad, no doubt about that.

Linux by default has compiling tools `make` and `gcc`; however in Docker we need to install both of them.

```bash
conda install -c conda-forge gcc=12.1.0
conda install -c anaconda make=4.2.1
```

In some cases, the c library `locale.h` needs to be re-linked.
Perl recognizes `xlocale.h`, which doesn't exist.
Here we are creating `xlocale.h` to link to the target `locale.h`.

```bash
cd /usr/include
sudo ln -s locale.h xlocale.h
```

## CPAN

Perl's package management is [CPAN](https://www.cpan.org/). Install the API from Conda:

```bash
conda install -c bioconda perl-app-cpanminus=1.7044
```

VEP requires Perl packages including `DBI` and `Try::Tiny`. Install them:

```bash
cpan DBI
cpan Try::Tiny
```

In addition, there is a MySQL package `DBD::mysql` required for VEP database connection.
`DBD::mysql` requires another dependency `mysql_config`, which was painful to install in Docker.
I ended up giving up, since we will be running VEP only in `--offline` mode without database connection.

## VEP

Install VEP in the folder `~/opt`

```bash
cd ~/opt
wget https://github.com/Ensembl/ensembl-vep/archive/release/106.zip
unzip 106.zip
rm 106.zip
```

You need to first `cd` into `ensembl-vep-release-106`,
otherwise some dependencies like `Bio` will not be installed properly within `ensembl-vep-release-106`.
This kind of hard-coded path of current directory was driving me crazy!

```bash
cd ensembl-vep-release-106
perl INSTALL.pl --NO_HTSLIB
```

You will see the warning message for `DBD::mysql` while installing VEP. Ignore it.

```bash
WARNING: DBD::mysql module not found. VEP can only run in offline (--offline) mode without DBD::mysql installed
```

There are several options during VEP installation. Here are my choices.

1. Do not cache the database. We will download it separately and run VEP in `--offline` mode.

```bash
The VEP can either connect to remote or local databases, or use local cache files.
Using local cache files is the fastest and most efficient way to run the VEP
Cache files will be stored in /home/linyc74/.vep
Do you want to install any cache files (y/n)? n
```

2. Do not download reference FASTA. We will use our own reference FASTA file.

```bash
The VEP can use FASTA files to retrieve sequence data for HGVS notations and reference sequence checks.
FASTA files will be stored in /home/linyc74/.vep
Do you want to install any FASTA files (y/n)? n
```

3. Install all plugins. These are source code, so they are small. Data resource files will be downloaded individually for each of them.

```bash
The VEP can use plugins to add functionality and data.
Plugins will be installed in /home/linyc74/.vep/Plugins
Do you want to install any plugins (y/n)? y
Plugins directory /home/linyc74/.vep/Plugins does not exists - do you want to create it (y/n)? y

The following plugins are available; which do you want (can specify multiple separated by spaces or 0 for all): 0
```

To reproduce the three settings automatically (in a Dockerfile),
we need to set `--AUTO` as `a` (API) and `p` (plugin),
and install all plugins by `--PLUGINS all`.
Also, include `--NO_UPDATE` so that the installation process
will not be disrupted by update check.
Disruption of installation will result in
missing perl modules (e.g. `Bio/EnsEMBL/Registry.pm`) and all plugins.

```bash
perl INSTALL.pl --AUTO ap --PLUGINS all --NO_HTSLIB --NO_UPDATE
```

## PATH

Make `vep` executable by adding to `PATH`.

```bash
export PATH=$PATH:$HOME/opt/ensembl-vep-release-106
```

In the end, there two folders related to VEP. If you want to cleanly uninstall VEP, just delete them.

- `~/opt/ensembl-vep-release-106`: Source code and build
- `~/.vep`: Config, cache, plugins... etc
