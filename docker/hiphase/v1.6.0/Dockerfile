FROM --platform=linux/amd64 ubuntu:22.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        wget curl git bzip2 ca-certificates tabix \
        libsm6 libxext6 libpcre3 libgl1 libpcre3-dev libssl-dev \
        zlib1g zlib1g-dev gfortran openssl \
        make build-essential cmake libncurses-dev \
        libbz2-dev liblzma-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda (x86_64)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_23.11.0-2-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Update conda
RUN conda update -n base -c defaults conda -y

# Install HiPhase
RUN conda install -c bioconda hiphase=1.6.0 -y

# Clean up to reduce Docker image size
RUN conda clean --all -y && \
    apt-get clean
