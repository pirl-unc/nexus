# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 ca-certificates \
    libsm6 libxext6 libpcre3 libgl1 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_24.1.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install HLAProfiler
RUN conda install -c bioconda hlaprofiler==1.0.5
RUN mkdir /opt/kraken/
RUN ln -s /opt/conda/share/kraken-ea-0.10.5ea.3-3 /opt/kraken

# Download HLAProfiler database
RUN wget https://github.com/ExpressionAnalysis/HLAProfiler/archive/refs/tags/v1.0.0-db_only.tar.gz && \
    mkdir /opt/HLAProfiler/ && \
    mv v1.0.0-db_only.tar.gz /opt/HLAProfiler/v1.0.0-db_only.tar.gz && \
    cd /opt/HLAProfiler/ && \
    tar -xvf v1.0.0-db_only.tar.gz && \
    cd /opt/HLAProfiler/HLAProfiler-1.0.0-db_only/hla_database && \
    wget https://github.com/ExpressionAnalysis/HLAProfiler/releases/download/v1.0.0-db_only/database.idx && \
    wget https://github.com/ExpressionAnalysis/HLAProfiler/releases/download/v1.0.0-db_only/database.kdb

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean