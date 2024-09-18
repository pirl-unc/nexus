# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

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
    tzdata  \
    wget \
    zlib1g \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py312_24.7.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update conda

# Install R
RUN conda config --set solver classic
RUN conda install -c conda-forge r-base==4.4.1
RUN conda install -c conda-forge libcurl==8.10.0

# Install Bambu
RUN Rscript -e "install.packages('optparse', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('BiocManager', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "BiocManager::install('rtracklayer')"
RUN Rscript -e "BiocManager::install('BSgenome')"
RUN Rscript -e "BiocManager::install('GenomicFeatures')"
RUN Rscript -e "BiocManager::install('bambu')"
ADD run_bambu.R /opt/bambu/run_bambu.R
RUN chmod +x /opt/bambu/run_bambu.R

# Clean up to reduce Docker image size
RUN apt-get clean
