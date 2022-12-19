FROM atgenomix.azurecr.io/0_runtime/seqslab_runtime-1.5_ubuntu-22.04_java11

#### Basic image utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    unzip \
    less \
    bwa=0.7.17-6 \
    autoconf \
    automake \
    gcc \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    git \
    unzip \
    cmake

#remove old tabix version, new tabix will be installed when install htslib v1.11
RUN apt-get remove -y tabix

# Install pandas
RUN python3 -m pip install pandas==1.3.5

# Link python3 to python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt install --no-install-recommends r-base=4.2.2.20221110-1.2204.0 -y

# Install FastQC v0.11.9
RUN cd /usr/local/seqslab \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && rm fastqc_v0.11.9.zip \
    && cd FastQC/ \
    && chmod 755 fastqc \
       ln -s `pwd`/fastqc /usr/local/bin/fastqc

# Install cutadapt v4.2
RUN python3 -m pip install cutadapt==4.2

# Install Trim_Galore v0.6.6
RUN cd /usr/local/seqslab \
    && curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz \
    && tar xvzf trim_galore.tar.gz \
    && rm trim_galore.tar.gz \
    && cd TrimGalore-0.6.6 \
    && ln -s `pwd`/trim_galore /usr/local/bin/trim_galore

# Install htslib v1.11, samtools v1.11, bcftools v1.8
RUN cd /usr/local/seqslab \
    && wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 \
    && wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
    && wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2 \
    && tar jxvf htslib-1.11.tar.bz2 \
    && tar jxvf samtools-1.11.tar.bz2 \
    && tar jxvf bcftools-1.8.tar.bz2 \
    && rm htslib-1.11.tar.bz2 samtools-1.11.tar.bz2 bcftools-1.8.tar.bz2 \
    && cd htslib-1.11 \
    && ./configure && make && make install && cd .. && rm -r htslib-1.11 \
    && cd bcftools-1.8 \
    && ./configure && make && make install && cd .. && rm -r bcftools-1.8 \
    && cd samtools-1.11 \
    && ./configure && make && make install && cd .. && rm -r samtools-1.11

# Install gatk-4.2.4.1
RUN cd /usr/local/seqslab \
    && wget https://github.com/broadinstitute/gatk/releases/download/4.2.4.1/gatk-4.2.4.1.zip \
    && unzip gatk-4.2.4.1.zip \
    && rm gatk-4.2.4.1.zip \
    && cd gatk-4.2.4.1 \
    && ln -s `pwd`/gatk /usr/local/bin/gatk

# Install MuSE
RUN cd /usr/local/seqslab \
    && wget https://github.com/wwylab/MuSE/archive/refs/tags/v1.0-rc.tar.gz \
    && tar zxvf v1.0-rc.tar.gz \
    && rm v1.0-rc.tar.gz \
    && cd MuSE-1.0-rc \
    && make \
    ln -s `pwd`/MuSE /usr/local/bin/muse

# Install VarDict
RUN apt-get install -y dos2unix \
    && cd /usr/local/seqslab \
    && git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git \
    && cd VarDictJava \
    && dos2unix gradlew \
    && ./gradlew distTar \
    && tar -xf build/distributions/VarDict-1.8.3.tar \
    && cp -r VarDict-1.8.3 .. \
    && cd .. \
    && rm -rf VarDictJava \
    && dos2unix ./VarDict-1.8.3/bin/*.pl \
    && dos2unix ./VarDict-1.8.3/bin/*.R \
    && ln -sr ./VarDict-1.8.3/bin/* /usr/local/bin

# Install LoFreq
RUN cd /usr/local/seqslab \
    && wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz \
    && tar -xzf lofreq_star-2.1.5_linux-x86-64.tgz \
    && rm lofreq_star-2.1.5_linux-x86-64.tgz \
    && mv lofreq_star-2.1.5_linux-x86-64/bin/* /usr/local/bin

# Install somatic-sniper
RUN cd /usr/local/seqslab \
    && wget https://github.com/genome/somatic-sniper/archive/refs/tags/v1.0.5.0.tar.gz \
    && tar zxvf v1.0.5.0.tar.gz \
    && rm v1.0.5.0.tar.gz \
    && cd somatic-sniper-1.0.5.0 \
    && mkdir somatic-sniper/build \
    && cd somatic-sniper/build \
    && cmake ../ \
    && make deps \
    && make -j \
    && ln -s `pwd`/bin/bam-somaticsniper /usr/local/bin

# Install bedtools
RUN cd /usr/local/seqslab \
    && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz \
    && tar -zxvf bedtools-2.30.0.tar.gz \
    && rm bedtools-2.30.0.tar.gz \
    && cd bedtools2 \
    && make \
    && ln -sr ./bin/* /usr/local/bin

# Install VarScan
COPY env/VarScan.v2.3.7.jar /usr/local/seqslab
RUN cd /usr/local/seqlab \
    && java -jar VarScan.v2.3.7.jar

# Install SnpEff
#RUN cd ~ \
#    && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
#    && unzip snpEff_latest_core.zip \
#    && rm snpEff_latest_core.zip \
#ENV PATH=$PATH:/snpEff/exec

# Install VEP
#RUN cd /usr/local/seqslab \
#    && wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/106.0.tar.gz \
#    && tar -zxvf 106.0.tar.gz \
#    && rm 106.0.tar.gz \
#    && cd ensembl-vep-release-106.0 \