FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive

######### dependencies
RUN apt-get -y update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
    less \
    bc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install MetONTIIME
WORKDIR /home/

RUN wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2023.9-py38-linux-conda.yml
RUN conda env create -n MetONTIIME_env --file qiime2-amplicon-2023.9-py38-linux-conda.yml
RUN rm qiime2-amplicon-2023.9-py38-linux-conda.yml

RUN conda install -n MetONTIIME_env -c bioconda seqtk
RUN /opt/conda/envs/MetONTIIME_env/bin/python -m pip install nanofilt

ENV PATH /opt/conda/envs/MetONTIIME_env/bin:$PATH
ENV CONDA_PREFIX /opt/conda/envs/MetONTIIME_env
