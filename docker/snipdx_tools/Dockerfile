#################################################################
# Dockerfile for snipdx pipeline (ver 2.1)
#
# Ubuntu             20.04
# R	                 3.6.1
# Python               3.8
# Facets            0.5.14
# pctGCdata          0.2.0
# reshape2           0.3.0 (version not specified when pulling from repo)
# snp-pileup        0.5.14
# Python packages:
# statistics, Pandas, PySam, numpy, matplotlib, seaborn
#################################################################

# Set the base image to Ubuntu.
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# Update and get basic dependencies
RUN apt-get -y update && apt-get install -y libxml2-dev  libcurl4-gnutls-dev libssl-dev texlive-base texinfo wget python3.8

# get miniconda locally to pull the rest of the dependencies in a simple and more reliable way than listing them
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh \
      && mkdir /root/.conda \
      && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b \
      &&  rm -f Miniconda3-py38_4.10.3-Linux-x86_64.sh

# install snp-pileup via conda through the bioconda channel
RUN conda install -c bioconda snp-pileup=0.5.14=ha52163a_2

# install r-facets via conda through the bioconda channel
RUN conda install -c bioconda r-facets=0.5.14=r36he991be0_2

# install R package RcppRoll
RUN R -e "install.packages(c('RcppRoll', 'optparse', 'reshape2'), dependencies=TRUE,  repos='http://cran.rstudio.com/')"

# install python packages
RUN pip install statistics pysam pandas

# additional packages for chromosome plots
RUN R -e "install.packages(c('tidyverse'), dependencies=TRUE,  repos='http://cran.rstudio.com/')"
RUN pip install numpy matplotlib seaborn
