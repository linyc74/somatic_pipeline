FROM continuumio/miniconda3:4.10.3

RUN conda install -c bioconda \
    trim-galore=0.6.6 \
    bwa=0.7.17 \
    samtools=1.11 \
    gatk4=4.2.4.1 \
 && conda install -c conda-forge \
    unzip=6.0 \
 && conda clean --all --yes

RUN pip install --no-cache-dir \
    pandas==1.3.4

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
 && unzip snpEff_latest_core.zip \
 && rm snpEff_latest_core.zip \
 && snpeff download GRCh38.99

ENV PATH /snpEff/exec:$PATH

COPY ./gatk_pipeline/* /gatk_pipeline/gatk_pipeline/
COPY ./__main__.py /gatk_pipeline/
WORKDIR /
