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
    tzdata \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py312_24.11.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
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

# Install minimap2
RUN conda install -c conda-forge libgcc-ng==14.2.0
RUN conda install -c bioconda minimap2==2.28

# Install abPOA
RUN wget https://github.com/yangao07/abPOA/releases/download/v1.4.1/abPOA-v1.4.1.tar.gz && \
    tar -xzf abPOA-v1.4.1.tar.gz && \
    cd abPOA-v1.4.1/ && \
    make

# Add abPOA to PATH
ENV PATH /abPOA-v1.4.1/bin:$PATH

# Install Mandarlorion 4.5
RUN conda install -c conda-forge mamba
RUN mamba install numpy Cython
RUN mamba install -c bioconda mappy
RUN wget https://github.com/christopher-vollmers/Mandalorion/archive/refs/tags/v4.5.tar.gz && \
    tar -xf v4.5.tar.gz

# Add Mandarlorion to PATH
ENV PATH /Mandalorion-4.5:$PATH

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean
