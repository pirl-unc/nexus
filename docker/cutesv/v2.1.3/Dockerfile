FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# TARGETARCH is set automatically by buildx: "amd64" or "arm64"
ARG TARGETARCH

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        wget curl git bzip2 ca-certificates \
        libpcre3 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
        make build-essential cmake libncurses-dev \
        libbz2-dev liblzma-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda — pick installer matching the target arch
RUN case "${TARGETARCH}" in \
        amd64) MINICONDA_ARCH="x86_64" ;; \
        arm64) MINICONDA_ARCH="aarch64" ;; \
        *) echo "Unsupported arch: ${TARGETARCH}" && exit 1 ;; \
    esac && \
    wget --quiet "https://repo.anaconda.com/miniconda/Miniconda3-py39_24.3.0-0-Linux-${MINICONDA_ARCH}.sh" -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Update conda
RUN conda update -n base -c defaults conda -y

# Install cuteSV and runtime deps
RUN conda install -c conda-forge libgcc-ng==13.2.0 -y && \
    conda install -c bioconda cutesv==2.1.3 -y

# Clean up to reduce image size
RUN conda clean --all -y && \
    apt-get clean