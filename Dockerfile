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

# install cnvkit, the pip used is in 'somatic' env
RUN conda install -c anaconda -n somatic pomegranate=0.14.4 \
 && conda clean --all --yes \
 && pip install --upgrade pip \
 && pip install --no-cache-dir cnvkit==0.9.9

# copy source code
COPY somatic_pipeline/* /somatic_pipeline/somatic_pipeline/
COPY ./__main__.py /somatic_pipeline/
WORKDIR /
