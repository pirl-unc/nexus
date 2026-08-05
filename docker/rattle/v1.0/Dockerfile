FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    bzip2 \
    ca-certificates \
    gzip \
    libpcre3 \
    libpcre3-dev \
    libssl-dev \
    zlib1g \
    zlib1g-dev \
    gfortran \
    openssl \
    make \
    build-essential \
    cmake \
    libncurses-dev \
    ncurses-dev \
    libbz2-dev \
    lzma-dev \
    liblzma-dev \
    libssl-dev \
    autoconf \
    automake \
    libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_24.3.0-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Instlal RATTLE
RUN conda install -c bioconda rattle==1.0

# Install pandas and pysam
RUN conda install -c anaconda pandas==2.2.3
RUN pip install pysam==0.22.1

ADD run_rattle_clustered.py /run_rattle_clustered.py
RUN chmod +x /run_rattle_clustered.py

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean