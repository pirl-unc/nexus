# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages
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

# Install openjdk
RUN conda install conda-forge::openjdk==22.0.1

# Download picard.jar
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar && \
    mkdir -p /opt/picard/ && \
    mv picard.jar /opt/picard/picard.jar && \
    chmod +x /opt/picard/picard.jar

# Add samtools to PATH
ENV PATH /opt/picard:$PATH

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean