FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive
ENV SILVA_VERSION=132
ENV GRGENES_VERSION=13_8

########## generate working directories
RUN mkdir /home/tools \
  && mkdir /home/references


######### dependencies
RUN apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install MetONTIIME
WORKDIR /home/tools/

RUN git clone https://github.com/MaestSi/MetONTIIME.git
WORKDIR /home/tools/MetONTIIME
RUN chmod 755 * \
  && ./install.sh

RUN sed -i 's/PIPELINE_DIR <- .*/PIPELINE_DIR <- \"\/home\/tools\/MetONTIIME\/\"/' config_MinION_mobile_lab.R
RUN sed -i 's/MINICONDA_DIR <- .*/MINICONDA_DIR <- \"\/opt\/conda\/\"/' config_MinION_mobile_lab.R


############################################################## download references
########################### SILVA
RUN mkdir /home/references/SILVA_${SILVA_VERSION}
WORKDIR /home/references/SILVA_${SILVA_VERSION}

RUN apt-get update -qq \
    && apt-get install --no-install-recommends -y \
    unzip \
    nano \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_${SILVA_VERSION}_release.zip \
  && unzip /home/references/SILVA_${SILVA_VERSION}/Silva_${SILVA_VERSION}_release.zip \
  && mv /home/references/SILVA_${SILVA_VERSION}/SILVA_${SILVA_VERSION}_QIIME_release/rep_set/rep_set_16S_only/99/silva_${SILVA_VERSION}_99_16S.fna /home/references/SILVA_${SILVA_VERSION}/silva_${SILVA_VERSION}_99_16S.fna \
  && mv /home/references/SILVA_${SILVA_VERSION}/SILVA_${SILVA_VERSION}_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt /home/references/SILVA_${SILVA_VERSION}/taxonomy_7_levels.txt \
  && rm -rf /home/references/SILVA_${SILVA_VERSION}/SILVA_${SILVA_VERSION}_QIIME_release \
  && rm -rf /home/references/SILVA_${SILVA_VERSION}/Silva_${SILVA_VERSION}_release.zip \
  && rm -rf /home/references/SILVA_${SILVA_VERSION}/__MACOSX

########################### Green Genes
RUN mkdir /home/references/GreenGenes
WORKDIR /home/references/GreenGenes

RUN wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_${GRGENES_VERSION}_otus.tar.gz \
  && tar xvzf /home/references/GreenGenes/gg_${GRGENES_VERSION}_otus.tar.gz \
  && mv /home/references/GreenGenes/gg_${GRGENES_VERSION}_otus/rep_set/99_otus.fasta /home/references/GreenGenes/99_otus.fasta \
  && mv /home/references/GreenGenes/gg_${GRGENES_VERSION}_otus/taxonomy/99_otu_taxonomy.txt /home/references/GreenGenes/99_otu_taxonomy.txt \
  && rm -rf /home/references/GreenGenes/gg_${GRGENES_VERSION}_otus \
  && rm -rf /home/references/GreenGenes/gg_${GRGENES_VERSION}_otus.tar.gz

WORKDIR /home/

