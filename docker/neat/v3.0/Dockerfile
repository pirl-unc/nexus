# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

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

# Download and Install Miniconda (Python 3.7)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add Miniconda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda -y

# Install NEAT dependencies (pin versions, add -y)
RUN conda install -c bioconda seqkit==2.12.0
RUN conda install -y numpy==1.19.5 && \
    conda install -y -c conda-forge biopython==1.78 && \
    conda install -y matplotlib==3.3.4 && \
    conda install -y pandas==1.2.1 && \
    pip install pysam==0.16.0.1 && \
    pip install matplotlib-venn==0.11.6

# Install NEAT
RUN mkdir /opt/neat/
RUN cd /opt/neat/ && \
    wget https://github.com/zstephens/neat-genreads/archive/refs/tags/3.0.tar.gz && \
    tar --strip-components=1 -xvf 3.0.tar.gz
RUN chmod +x /opt/neat/gen_reads.py
RUN chmod -R +x /opt/neat/source/*

ENV PATH="/opt/neat:${PATH}"

# Clean up to reduce Docker image size
RUN conda clean -a -y && \
    apt-get clean
