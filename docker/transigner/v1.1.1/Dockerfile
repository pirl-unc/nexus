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
    tzdata  \
    wget \
    zlib1g \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_24.7.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update conda

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
RUN conda install -c bioconda minimap2==2.22

# Install transigner
RUN conda config --set solver classic
RUN conda install pandas==2.2.2
RUN pip install transigner==1.1.1

# Clean up to reduce Docker image size
RUN apt-get clean
