FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 ca-certificates \
    libsm6 libxext6 libpcre3 libgl1 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev tabix \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_23.11.0-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Install isoLASER's pinned dependencies first so pip treats them as satisfied
# and does not change them when isoLASER is installed.
RUN pip install --no-cache-dir \
        multiprocess \
        numpy==1.24.4 \
        scipy==1.10.1 \
        pandas==2.0.3 \
        pysam==0.16.0.1 \
        pytabix==0.1 \
        pyfaidx==0.5.8 \
        biopython==1.76 \
        vcfpy==0.13.8 \
        HTSeq==0.12.4 \
        networkx==3.1 \
        scikit-learn==1.3.2

RUN pip install --no-cache-dir isolaser==0.0.3

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
RUN conda install -c conda-forge libgcc-ng==14.2.0
RUN conda install -c bioconda minimap2==2.28

# Clean up to reduce Docker image size
RUN conda clean --all --yes
RUN apt-get clean