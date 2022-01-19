FROM continuumio/miniconda3:4.10.3

RUN conda create -n gatk \
 && conda install -c bioconda -n gatk \
    trim-galore=0.6.6 \
    bwa=0.7.17 \
    samtools=1.11 \
    gatk4=4.2.4.1 \
 && conda install -c anaconda -n gatk \
    pandas=1.3.5 \
 && conda clean --all --yes

ENV PATH /opt/conda/envs/gatk/bin:$PATH

RUN conda install -c conda-forge unzip=6.0 \
 && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
 && unzip snpEff_latest_core.zip \
 && rm snpEff_latest_core.zip

ENV PATH /snpEff/exec:$PATH

RUN snpeff download -verbose GRCh38.99

COPY ./gatk_pipeline/* /gatk_pipeline/gatk_pipeline/
COPY ./__main__.py /gatk_pipeline/
WORKDIR /
