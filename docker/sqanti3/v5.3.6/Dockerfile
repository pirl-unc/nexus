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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_24.7.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Update Conda
RUN conda update -y conda

# Install SQANTI3 5.3.6
RUN wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.3.6/SQANTI3_v5.3.6.zip && \
    mkdir sqanti3 && \
    unzip SQANTI3_v5.3.6.zip -d sqanti3/

# Add SQANTI3 to PATH
ENV PATH /sqanti3:$PATH

# Install SQANTI3 environment
RUN conda env update --name base --file sqanti3/SQANTI3.conda_env.yml

# Install Minimap2 v2.28
RUN conda install -c bioconda minimap2==2.28

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean