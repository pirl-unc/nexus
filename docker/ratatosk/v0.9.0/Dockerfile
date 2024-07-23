FROM ubuntu:18.04
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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_24.1.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install necessary packages
RUN conda install -c conda-forge libgcc-ng==12.2.0

RUN wget https://github.com/DecodeGenetics/Ratatosk/archive/refs/tags/v0.9.0.tar.gz && \
    tar -xvf v0.9.0.tar.gz
RUN wget https://github.com/pmelsted/bifrost/archive/refs/tags/v1.2.1.tar.gz && \
    tar -xvf v1.2.1.tar.gz
RUN mv bifrost-1.2.1/* Ratatosk-0.9.0/Bifrost/
RUN cd Ratatosk-0.9.0/ && \
    mkdir build &&  \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    chmod +x /Ratatosk-0.9.0/build/src/Ratatosk

# Add Ratatosk to PATH
ENV PATH=/Ratatosk-0.9.0/build/src/:$PATH

# Clean up to reduce Docker image size
RUN conda clean -all