# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 ca-certificates  \
    libpcre3 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xvjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18/ && \
    mkdir bin && \
    ./configure --prefix=$(pwd) && \
    make && \
    make install

# Add samtools to PATH
ENV PATH /samtools-1.18:$PATH

# Install liqa
RUN pip install liqa==1.3.4

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean