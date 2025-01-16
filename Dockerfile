FROM continuumio/miniconda3:24.11.1-0

RUN apt-get update

RUN conda create -n somatic \
 && conda install -n somatic -c anaconda \
    pandas=2.1.4 \
 && conda install -n somatic -c conda-forge \
    tbb=2020.2 \
    unzip=6.0

RUN conda install -n somatic -c bioconda \
    trim-galore=0.6.6 \
    bwa=0.7.18 \
    samtools=1.21 \
    gatk4=4.6.1.0 \
    bowtie2=2.5.4 \
    muse=1.0 \
    varscan=2.3.7 \
    bcftools=1.8 \
    vcf2maf=1.6.21

# Activate somatic env but do not prioritize over the system env (for the perl used by VEP)
ENV PATH=$PATH:/opt/conda/envs/somatic/bin

# System lib is prioritized over somatic env lib
ENV LD_LIBRARY_PATH=/lib/x86_64-linux-gnu:/opt/conda/envs/somatic/lib:$LD_LIBRARY_PATH

# Resolve bcftools dependency issue
ARG d=/opt/conda/envs/somatic/lib/
RUN ln -s ${d}libcrypto.so.1.1 ${d}libcrypto.so.1.0.0

# --- LoFreq ---
RUN wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz \
 && tar -xzf lofreq_star-2.1.5_linux-x86-64.tgz \
 && rm lofreq_star-2.1.5_linux-x86-64.tgz
ENV PATH=$PATH:/lofreq_star-2.1.5_linux-x86-64/bin

# --- SomaticSniper ---
RUN conda install -c bioconda -n somatic \
    bedtools=2.30.0 \
    somatic-sniper=1.0.5.0

# --- VarDict ---
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
ENV PATH=$PATH:/VarDict-1.8.3/bin

# --- VEP ---
# Perpare dependencies for VEP installation, a nightmare because every new build has different problems
RUN apt-get install -y \
    build-essential \
    perl \
    cpanminus \
    curl \
    libwww-perl \
 && cpanm DBI \
 && cpanm Try::Tiny
# Include "--NO_UPDATE" so that the installation process will not be disrupted by update check
# Disruption of installation will result in missing perl modules (e.g. Bio/EnsEMBL/Registry.pm) and plugins
RUN wget https://github.com/Ensembl/ensembl-vep/archive/release/106.zip \
 && unzip 106.zip \
 && rm 106.zip \
 && cd ensembl-vep-release-106 \
 && perl INSTALL.pl --AUTO ap --PLUGINS all --NO_HTSLIB --NO_UPDATE \
 && cd ..
ENV PATH=$PATH:/ensembl-vep-release-106

# --- PCGR ---
ENV PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v2.1.2/conda/env/lock/"
ENV PLATFORM="linux"
RUN conda create --name pcgr --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock
RUN conda create --name pcgrr --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock
ENV PATH=$PATH:/opt/conda/envs/pcgr/bin
RUN cp /opt/conda/envs/pcgr/bin/pcgrr.R /opt/conda/envs/somatic/bin/

RUN apt-get autoremove \
 && apt-get clean \
 && conda clean --all --yes

COPY somatic_pipeline/* /somatic_pipeline/somatic_pipeline/
COPY ./__main__.py /somatic_pipeline/
WORKDIR /
