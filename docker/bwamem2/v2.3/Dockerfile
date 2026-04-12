FROM ubuntu:20.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        wget curl git bzip2 ca-certificates \
        libsm6 libxext6 libgl1 \
        libpcre3 libpcre3-dev \
        libssl-dev openssl \
        zlib1g zlib1g-dev \
        gfortran \
        make build-essential cmake \
        libncurses-dev ncurses-dev \
        libbz2-dev lzma-dev liblzma-dev \
        autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py312_25.9.1-3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Install bwa-mem2
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda update -n base -c defaults conda -y && \
    conda install -c bioconda -c conda-forge bwa-mem2==2.3 samtools==1.21 -y && \
    conda clean -a -y

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean