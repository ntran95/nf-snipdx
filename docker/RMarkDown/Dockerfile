####################################################################
# Dockerfile for running RMarkDown scripts
#
# Ubuntu             20.04
# R			         4.1.0
####################################################################

# Set the base image to Ubuntu.
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive

# Update and get basic dependencies
RUN apt-get -y update && apt-get install -y libxml2-dev  libcurl4-gnutls-dev libssl-dev texlive-base \
    texinfo wget python3.8 gcc g++ pandoc gfortran fort77 gcc-multilib libreadline-dev \
    xorg-dev liblzma-dev libblas-dev fort77 libbz2-dev libcairo2-dev build-essential  \
    && apt-get clean

# install prec2
RUN wget https://github.com/PhilipHazel/pcre2/releases/download/pcre2-10.39/pcre2-10.39.tar.bz2 \
    && tar -xjvf pcre2-10.39.tar.bz2 \
    && cd pcre2-10.39 \
    && ./configure \
    && make \
    && make install \
    && cd ../ \
    && rm pcre2-10.39.tar.bz2

# install R 4.1
RUN wget -c https://cran.r-project.org/src/base/R-4/R-4.1.0.tar.gz \
    && tar -xf R-4.1.0.tar.gz \
    && cd R-4.1.0 \
    && ./configure --with-readline=no --with-x=no \
    && make -j9 \
    && make install \
    && cd ../ \
    && rm R-4.1.0.tar.gz

# install geneal R packages from CRAN
RUN R -e "install.packages(c('rmarkdown', 'readxl', 'openxlsx', 'reshape2', 'DT', 'tidyverse',  \
        'ggrepel', 'here', 'cowplot', 'purrr', 'data.table', 'binom', 'magrittr', 'BiocManager'), \
        dependencies=TRUE,  repos='http://cran.rstudio.com/')"

# install R packages from Bioconductor
RUN R -e "BiocManager::install(c('GenomicRanges'), dependencies=TRUE)"
