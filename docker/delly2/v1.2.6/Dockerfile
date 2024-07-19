# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Download Delly2
RUN mkdir /opt/delly2/
RUN wget https://github.com/dellytools/delly/releases/download/v1.2.6/delly_v1.2.6_linux_x86_64bit && \
    mv delly_v1.2.6_linux_x86_64bit /opt/delly2/delly && \
    chmod +x /opt/delly2/delly

# Add delly to PATH
ENV PATH /opt/delly2/:$PATH

# Install bcftools
RUN mkdir /opt/bcftools/
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xvjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/opt/bcftools/ && \
    make && \
    make install

# Add bcftools to PATH
ENV PATH /opt/bcftools/bin:$PATH

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

# Clean up to reduce Docker image size
RUN apt-get clean