FROM continuumio/miniconda3:4.10.3

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
    pandas=1.3.5 \
 && conda clean --all --yes

# for identical commands (e.g. pip), somatic overrides default environment
ENV PATH /opt/conda/envs/somatic/bin:$PATH

# bcftools dependency issue
ARG d=/opt/conda/envs/somatic/lib/
RUN ln -s ${d}libcrypto.so.1.1 ${d}libcrypto.so.1.0.0

# download and unzip snpeff
RUN conda install -c conda-forge unzip=6.0 \
 && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
 && unzip snpEff_latest_core.zip \
 && rm snpEff_latest_core.zip

# make snpeff executable
ENV PATH /snpEff/exec:$PATH

# download pre-build snpeff database
RUN snpeff download -verbose GRCh38.99

# dependency for cnvkit
ARG d=/opt/conda/envs/somatic/lib/
RUN conda install -c conda-forge -n somatic r-base=3.2.2 \
 && conda install -c bioconda -n somatic bioconductor-dnacopy=1.44.0 \
 && ln -s ${d}libreadline.so.8.1 ${d}libreadline.so.6 \
 && ln -s ${d}libncursesw.so.6.2 ${d}libncurses.so.5 \
 && conda install -c anaconda -n somatic pomegranate=0.14.4 \
 && conda clean --all --yes

# install cnvkit, the pip used is in 'somatic' env
RUN pip install --upgrade pip \
 && pip install --no-cache-dir cnvkit==0.9.9

# perl dependency for vep
# perl build must be "h470a237_0" to avoid bad version (hard-coded gcc path)
RUN conda install -c conda-forge -n somatic \
    perl=5.26.2=h470a237_0 \
    gcc=12.1.0 \
 && conda install -c anaconda -n somatic \
    make=4.2.1 \
 && conda install -c bioconda -n somatic \
    perl-app-cpanminus=1.7044 \
 && cpan DBI \
 && cpan Try::Tiny \
 && cpan Bio::EnsEMBL::Registry

# install vep
RUN wget https://github.com/Ensembl/ensembl-vep/archive/release/106.zip \
 && unzip 106.zip \
 && perl ensembl-vep-release-106/INSTALL.pl --AUTO p --PLUGINS all --NO_HTSLIB \
 && rm 106.zip

# make vep executable
ENV PATH /ensembl-vep-release-106:$PATH

# copy source code
COPY somatic_pipeline/* /somatic_pipeline/somatic_pipeline/
COPY ./__main__.py /somatic_pipeline/
WORKDIR /
