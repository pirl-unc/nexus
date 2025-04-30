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
    glibc-source \
    libglib2.0-0 \
    libc6 \
    libc6-dev \
    libc-dev \
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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install Talon
RUN pip install scanpy
RUN pip install pysam
RUN pip install cython
RUN conda install -c conda-forge "libcurl>=7.22.0,<8"
RUN conda install pandas
RUN conda install -c bioconda pyranges
RUN wget https://github.com/mortazavilab/TALON/archive/refs/tags/v6.0.1.tar.gz && \
    tar -xf v6.0.1.tar.gz && \
    cd TALON-6.0.1 && \
    pip install .

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean
