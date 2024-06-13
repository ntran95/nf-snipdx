FROM ubuntu:20.04

USER root

RUN apt-get -yq update
ENV DEBIAN_FRONTEND "noninteractive"

# Update and get basic dependencies
RUN apt-get -y update && apt-get install -y libxml2-dev  libcurl4-gnutls-dev libssl-dev texlive-base \
    texinfo wget python3.8 default-jre-headless build-essential python3-pip libgd-dev libbz2-dev liblzma-dev \
    && apt-get clean

RUN pip3 install pandas pysam
