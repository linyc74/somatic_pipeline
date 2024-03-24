FROM continuumio/miniconda3:4.12.0

# create conda env & install packages
RUN conda create -n somatic \
 && conda install -c conda-forge -n somatic \
    tbb=2020.2 \
 && conda install -c bioconda -n somatic \
    trim-galore=0.6.6 \
    bwa=0.7.17 \
    samtools=1.11 \
    gatk4=4.2.4.1 \
    bowtie2=2.3.5 \
    muse=1.0 \
    varscan=2.3.7 \
    bcftools=1.8 \
    vcf2maf=1.6.21 \
 && conda install -c anaconda -n somatic \
    pandas=1.3.5

# activate somatic env
# for identical commands (e.g. pip), somatic overrides default environment
ENV PATH /opt/conda/envs/somatic/bin:$PATH

# system lib is prioritized over somatic env lib
ENV LD_LIBRARY_PATH=/lib/x86_64-linux-gnu:/opt/conda/envs/somatic/lib:$LD_LIBRARY_PATH

# bcftools dependency issue
ARG d=/opt/conda/envs/somatic/lib/
RUN ln -s ${d}libcrypto.so.1.1 ${d}libcrypto.so.1.0.0

# --- lofreq ---
# download and untar lofreq
RUN wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz \
 && tar -xzf lofreq_star-2.1.5_linux-x86-64.tgz \
 && rm lofreq_star-2.1.5_linux-x86-64.tgz
# make lofreq executable
ENV PATH /lofreq_star-2.1.5_linux-x86-64/bin:$PATH

# --- somatic-sniper ---
RUN conda install -c bioconda -n somatic \
    bedtools=2.30.0 \
    somatic-sniper=1.0.5.0

# --- vardict ---
RUN apt-get update \
 && apt-get install -y dos2unix \
 && git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git \
 && cd VarDictJava \
 && dos2unix gradlew \
 && ./gradlew distTar \
 && tar -xf build/distributions/VarDict-1.8.3.tar \
 && cp -r VarDict-1.8.3 / \
 && cd .. \
 && rm -rf VarDictJava \
 && dos2unix /VarDict-1.8.3/bin/*.pl \
 && dos2unix /VarDict-1.8.3/bin/*.R
# make vardict executable
ENV PATH=/VarDict-1.8.3/bin:$PATH

# --- snpeff ---
# download and unzip snpeff
RUN conda install -c conda-forge unzip=6.0 \
 && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
 && unzip snpEff_latest_core.zip \
 && rm snpEff_latest_core.zip
# make snpeff executable
ENV PATH /snpEff/exec:$PATH

# --- vep ---
# perl dependency for vep
# perl build must be "5.26.2=h470a237_0" to avoid bad version (hard-coded gcc path)
#RUN conda install -c conda-forge -n somatic \
#    perl=5.26.2=h470a237_0 \
#    gcc=12.1.0 \
# && conda install -c anaconda -n somatic \
#    make=4.2.1 \
# && conda install -c bioconda -n somatic \
#    perl-app-cpanminus=1.7044 \
# && cpan DBI \
# && cpan Try::Tiny

RUN apt-get install -y \
    build-essential \
    perl \
    cpanminus \
 && cpanm DBI \
 && cpanm Try::Tiny

# install vep
# include "--NO_UPDATE" so that the installation process will not be disrupted by update check
# disruption of installation will result in missing perl modules (e.g. Bio/EnsEMBL/Registry.pm) and plugins 
RUN wget https://github.com/Ensembl/ensembl-vep/archive/release/106.zip \
 && unzip 106.zip \
 && rm 106.zip \
 && cd ensembl-vep-release-106 \
 && perl INSTALL.pl --AUTO ap --PLUGINS all --NO_HTSLIB --NO_UPDATE \
 && cd ..
# make vep executable
ENV PATH /ensembl-vep-release-106:$PATH

# clean up
RUN apt-get autoremove \
 && apt-get clean \
 && conda clean --all --yes

# copy source code
COPY somatic_pipeline/* /somatic_pipeline/somatic_pipeline/
COPY ./__main__.py /somatic_pipeline/
WORKDIR /
