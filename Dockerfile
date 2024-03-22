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



# --- snpeff ---
# download and unzip snpeff
RUN conda install -c conda-forge unzip=6.0 \
 && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
 && unzip snpEff_latest_core.zip \
 && rm snpEff_latest_core.zip

# make snpeff executable
ENV PATH /snpEff/exec:$PATH



# --- R ---
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
 && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
 && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
 && apt-get install -y --no-install-recommends \
    r-base-dev=4.0.4-1



# --- cnvkit ---
# the pip used is in 'somatic' env
RUN Rscript -e 'install.packages("BiocManager", version="3.16")' \
 && Rscript -e 'BiocManager::install("DNAcopy")' \
 && conda install -c anaconda -n somatic \
    pomegranate=0.14.4 \
 && pip install --upgrade pip \
 && pip install --no-cache-dir \
    cnvkit==0.9.9



# --- vep ---
# perl dependency for vep
# perl build must be "5.26.2=h470a237_0" to avoid bad version (hard-coded gcc path)
RUN conda install -c conda-forge -n somatic \
    perl=5.26.2=h470a237_0 \
    gcc=12.1.0 \
 && cp /usr/include/crypt.h /opt/conda/envs/somatic/lib/5.26.2/x86_64-linux-thread-multi/CORE/crypt.h \
 && conda install -c anaconda -n somatic \
    make=4.2.1 \
 && conda install -c bioconda -n somatic \
    perl-app-cpanminus=1.7044 \
 && cpan DBI \
 && cpan Try::Tiny \
 && cpan LWP::Simple

# install vep
# include "--NO_UPDATE" so that the installation process will not be disrupted by update check
# disruption of installation will result in missing perl modules (e.g. Bio/EnsEMBL/Registry.pm) and plugins 
RUN wget https://github.com/Ensembl/ensembl-vep/archive/release/110.zip \
 && unzip 110.zip \
 && rm 110.zip \
 && cd ensembl-vep-release-110 \
 && perl INSTALL.pl --AUTO ap --PLUGINS all --NO_HTSLIB --NO_UPDATE \
 && cd ..

# make vep executable
ENV PATH /ensembl-vep-release-110:$PATH



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
RUN apt-get install -y dos2unix \
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

ENV PATH=/VarDict-1.8.3/bin:$PATH



# clean up
RUN apt-get autoremove \
 && apt-get clean \
 && conda clean --all --yes



# copy source code
COPY somatic_pipeline/* /somatic_pipeline/somatic_pipeline/
COPY ./__main__.py /somatic_pipeline/
WORKDIR /
