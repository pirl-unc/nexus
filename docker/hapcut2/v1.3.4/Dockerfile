FROM --platform=linux/amd64 ubuntu:22.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        wget curl git bzip2 ca-certificates tabix \
        libpcre3 libpcre3-dev libssl-dev \
        zlib1g zlib1g-dev gfortran openssl \
        make build-essential cmake libncurses-dev \
        libbz2-dev liblzma-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda (x86_64)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_25.7.0-2-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Accept conda Terms of Service
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Configure channels per the bioconda recommendation (conda-forge highest priority).
# HapCUT2 1.3.4 depends on liblzma/libzlib that only live on conda-forge.
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Update conda
RUN conda update -n base conda -y

# Install HapCUT2
RUN conda install hapcut2==1.3.4 -y

# Clean up to reduce Docker image size
RUN conda clean --all -y && \
    apt-get clean
