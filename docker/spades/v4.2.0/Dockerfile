FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 time ca-certificates \
    libsm6 libxext6 libpcre3 libgl1 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_26.3.2-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

RUN mkdir /opt/spades/ && \
    cd /opt/spades/ && \
    wget https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz && \
    tar -xvf SPAdes-4.2.0-Linux.tar.gz

ENV PATH=/opt/spades/SPAdes-4.2.0-Linux/bin:$PATH

# Clean up to reduce Docker image size
RUN conda clean --all --yes
RUN apt-get clean