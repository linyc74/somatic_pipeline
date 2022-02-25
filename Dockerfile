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
 && conda install -c anaconda -n somatic \
    pandas=1.3.5 \
 && conda clean --all --yes

# make conda env available to environment PATH
ENV PATH /opt/conda/envs/somatic/bin:$PATH

# download and unzip snpeff
RUN conda install -c conda-forge unzip=6.0 \
 && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
 && unzip snpEff_latest_core.zip \
 && rm snpEff_latest_core.zip

# make snpeff executable
ENV PATH /snpEff/exec:$PATH

# download pre-build snpeff database
RUN snpeff download -verbose GRCh38.99

# bcftools dependency issue
ARG d=/opt/conda/envs/somatic/lib/
RUN ln -s ${d}libcrypto.so.1.1 ${d}libcrypto.so.1.0.0

# copy source code
COPY somatic_pipeline/* /somatic_pipeline/somatic_pipeline/
COPY ./__main__.py /somatic_pipeline/
WORKDIR /
