# Set the base image to Ubuntu.
FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && apt-get install -y wget unzip default-jre && apt-get clean

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

RUN unzip snpEff_latest_core.zip
