FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 unzip ca-certificates  \
    libpcre3 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_24.1.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

RUN conda update -n base -c defaults conda

# Download and unpack HLAminer v1.4.2
WORKDIR /opt
RUN wget --quiet https://github.com/BirolLab/HLAminer/releases/download/v1.4.2/hlaminer-1.4.2.tar.gz && \
    tar -xvf hlaminer-1.4.2.tar.gz && \
    rm hlaminer-1.4.2.tar.gz
RUN ln -s /opt/hlaminer-1.4.2/HLAminer-1.4 /opt/hlaminer

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xvjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18/ && \
    mkdir bin && \
    ./configure --prefix=$(pwd) && \
    make && \
    make install

# Add samtools to PATH
# `make install --prefix=/opt/samtools-1.18` (from inside that dir) installs
# the samtools binary to /opt/samtools-1.18/bin/samtools.
ENV PATH /opt/samtools-1.18/bin:$PATH

# Install minimap2 for long-read alignment against HLAminer's HLA reference
# (used by the hla_typing_hlaminer-lr-rna subworkflow).
RUN conda install -c bioconda minimap2=2.28 -y

RUN conda clean -a -y && \
    apt-get clean