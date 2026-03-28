# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Install packages
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
      wget \
      curl \
      git \
      openssh-client \
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
      pkg-config \
      libncurses-dev \
      ncurses-dev \
      libbz2-dev \
      lzma-dev \
      liblzma-dev \
      autoconf \
      automake \
      libcurl4-openssl-dev \
      liblzma-dev \
      libdeflate-dev \
      libdeflate0 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py312_24.11.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install LongTR
RUN mkdir /opt/longtr/ && \
    cd /opt/longtr/ && \
    wget https://github.com/gymrek-lab/LongTR/archive/refs/tags/v1.2.tar.gz && \
    tar -xf v1.2.tar.gz && \
    cd LongTR-1.2/ && \
    git config --global url."https://github.com/".insteadOf git@github.com: && \
    sed -i 's/-lcrypto/-lcrypto -ldeflate/' Makefile && \
    make

ENV PATH=/opt/longtr/LongTR-1.2:$PATH

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean
