# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    build-essential \
    bzip2 \
    ca-certificates \
    cmake \
    curl \
    gcc-12 \
    gfortran \
    git \
    libarchive-dev \
    libncurses-dev \
    libpng-dev \
    libpcre3 \
    libpcre3-dev \
    libssl-dev \
    libbz2-dev \
    lzma-dev \
    liblzma-dev \
    libssl-dev \
    libpcre2-dev \
    libcurl4-openssl-dev \
    make \
    ncurses-dev \
    openssl \
    openmpi-common \
    libopenmpi-dev \
    samtools \
    tabix \
    tzdata  \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update conda

# Install REDItools2
RUN wget https://github.com/tflati/reditools2.0/archive/master.zip && \
    mkdir /opt/reditools2/ &&  \
    mv master.zip /opt/reditools2/ && \
    cd /opt/reditools2/ &&  \
    unzip master.zip && \
    cd reditools2.0-master/ && \
    bash install.sh
RUN sed -i '1i#!/usr/bin/env python' /opt/reditools2/reditools2.0-master/src/cineca/reditools_table_to_bed.py
RUN sed -i '1i#!/usr/bin/env python' /opt/reditools2/reditools2.0-master/src/cineca/annotate_with_DNA.py
RUN chmod +x /opt/reditools2/reditools2.0-master/src/cineca/reditools_table_to_bed.py && \
    chmod +x /opt/reditools2/reditools2.0-master/src/cineca/annotate_with_DNA.py && \
    ln -s /opt/reditools2/reditools2.0-master/src/cineca/reditools_table_to_bed.py /usr/local/bin/reditools_table_to_bed.py && \
    ln -s /opt/reditools2/reditools2.0-master/src/cineca/annotate_with_DNA.py /usr/local/bin/annotate_with_DNA.py

# Add REDItools2 to PATH
ENV PATH /opt/reditools2/reditools2.0-master/src/cineca/:$PATH

# Install REDItools1
RUN wget https://github.com/BioinfoUNIBA/REDItools/archive/refs/tags/REDItoolsFamily_V1.0.tar.gz && \
    mkdir /opt/reditools/ &&  \
    mv REDItoolsFamily_V1.0.tar.gz /opt/reditools/ && \
    cd /opt/reditools/ &&  \
    tar -xvf REDItoolsFamily_V1.0.tar.gz && \
    cd REDItools-REDItoolsFamily_V1.0/ && \
    python setup.py install
RUN pip install pysam==0.20.0

# Clean up to reduce Docker image size
RUN apt-get clean
