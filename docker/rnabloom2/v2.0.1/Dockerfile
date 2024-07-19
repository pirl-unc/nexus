# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 ca-certificates  \
    libpcre3 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl openjdk-11-jdk openjdk-11-jre \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_24.3.0-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xvjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18/ && \
    mkdir bin && \
    ./configure --prefix=$(pwd) && \
    make && \
    make install

# Add samtools to PATH
ENV PATH /samtools-1.18:$PATH

# Install minimap2
RUN conda install -c bioconda minimap2==2.22

# Install racon
RUN conda install -c conda-forge libzlib==1.2.11
RUN conda install -c bioconda racon==1.5.0

# Install ntcard
RUN mkdir /opt/ntcard/
RUN wget https://github.com/bcgsc/ntCard/releases/download/1.2.2/ntcard-1.2.2.tar.gz && \
    tar -xvf ntcard-1.2.2.tar.gz && \
    cd ntcard-1.2.2/ && \
    ./configure --prefix=/opt/ntcard/ && \
    make && \
    make install

# Add ntcard to PATH
ENV PATH /opt/ntcard/bin:$PATH

# Download rna-bloom2
RUN mkdir /opt/rnabloom2/ && \
    cd /opt/rnabloom2/ && \
    wget https://github.com/bcgsc/RNA-Bloom/releases/download/v2.0.1/rnabloom_v2.0.1.tar.gz && \
    tar -xvf rnabloom_v2.0.1.tar.gz && \
    mv RNA-Bloom_v2.0.1/RNA-Bloom.jar RNA-Bloom.jar && \
    chmod +x /opt/rnabloom2/RNA-Bloom.jar

# Add RNA-Bloom.jar to PATH
ENV PATH /opt/rnabloom2/:$PATH

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean