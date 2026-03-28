# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py312_24.11.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install hifiasm
RUN mkdir -p /opt/hifiasm/ && \
    cd /opt/hifiasm/ && \
    wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.25.0.tar.gz && \
    tar -xvf 0.25.0.tar.gz && \
    cd /opt/hifiasm/hifiasm-0.25.0/ && \
    make && \
    cp /opt/hifiasm/hifiasm-0.25.0/hifiasm /usr/local/bin/ && \
    chmod 755 /usr/local/bin/hifiasm

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean
