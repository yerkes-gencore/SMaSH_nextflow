# syntax=docker/dockerfile:1

## NOTE: This is a work in progress and is definitely not functional yet!


# Set the base image to python
FROM python:3.10-buster

# numpy==1.22.1
# pysam==0.18.0
# scipy==1.7.3


# Install STAR
ARG STAR_VERSION=2.7.9a

ENV STAR_DEPENDS gcc g++ make wget zlib1g-dev unzip

RUN set -ex

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${STAR_DEPENDS} && \
    apt-get clean && \
    g++ --version && \
    cd /home && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip && \
    unzip ${STAR_VERSION}.zip && \
    cd STAR-${STAR_VERSION}/source && \
    make STARstatic && \
    mkdir /home/bin && \
    cp STAR /home/bin && \
    cd /home && \
    'rm' -rf STAR-${STAR_VERSION} && \
    apt-get --purge autoremove -y ${STAR_DEPENDS}

ENV PATH /home/bin:${PATH}

# Install samtools
ARG STAR_VERSION=2.7.9a

# ENV SAMTOOLS_DEPENDS 

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    g++ --version && \
    cd /home && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip && \
    unzip ${STAR_VERSION}.zip && \
    cd STAR-${STAR_VERSION}/source && \
    make STARstatic && \
    mkdir /home/bin && \
    cp STAR /home/bin && \
    cd /home && \
    'rm' -rf STAR-${STAR_VERSION} && \
    apt-get --purge autoremove -y ${PACKAGES}

ENV PATH /home/bin:${PATH}

# Install samtools
RUN apt-get update & apt-get install --yes --no-install-recommends \
    locales \
    vim-tiny \
    python


