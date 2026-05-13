FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

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
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses-dev \
    libpcre2-dev \
    libpcre3 \
    libpcre3-dev \
    libpng-dev \
    libssl-dev \
    make \
    openssl \
    tzdata \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_24.11.1-0-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

RUN conda update -n base -c defaults conda -y

RUN pip install --no-cache-dir pysam==0.23.0 && \
    conda install -y -c conda-forge rust==1.88.0 polars==1.26.0 pyarrow==19.0.1 && \
    conda install -y -c anaconda pandas==2.2.3 && \
    conda clean --all --yes

RUN mkdir -p /opt/exacto/ && \
    cd /opt/exacto/ && \
    wget https://github.com/pirl-unc/exacto/releases/download/v0.4.6a1/exacto-0.4.6a1.tar.gz && \
    pip install --no-cache-dir exacto-0.4.6a1.tar.gz && \
    conda clean --all --yes