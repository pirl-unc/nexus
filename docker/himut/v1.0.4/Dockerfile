# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    bzip2 \
    ca-certificates \
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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.9.1-3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Install himut
RUN pip install biopython==1.78
RUN pip install cyvcf2==0.30.18
RUN pip install matplotlib==3.5.0
RUN pip install mizani==0.8.1
RUN pip install natsort==8.0.0
RUN pip install numpy==1.20.2
RUN pip install plotnine==0.9.0
RUN pip install psutil==5.8.0
RUN pip install pyfastx==0.8.4
RUN pip install pysam==0.18.0
RUN pip install pytabix==0.1
RUN pip install scipy==1.9.0
RUN wget https://github.com/sjin09/himut/archive/refs/tags/v1.0.4.tar.gz && \
    pip install v1.0.4.tar.gz

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean