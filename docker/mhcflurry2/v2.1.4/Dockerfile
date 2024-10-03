# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

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
    gcc-12 \
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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_24.7.1-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update conda

# Install mhcflurry
RUN conda config --set solver classic
RUN conda install conda-forge::six==1.16.0
RUN conda install pandas==2.2.2
RUN conda install conda-forge::tensorflow==2.17.0
RUN conda install conda-forge::keras==3.5.0
RUN conda install conda-forge::appdirs==1.4.4
RUN conda install anaconda::scikit-learn==1.5.1
RUN conda install anaconda::pyyaml==6.0.1
RUN conda install conda-forge::tqdm==4.66.5
RUN pip install mhcgnomes==1.8.6
RUN pip install mhcflurry
RUN mhcflurry-downloads fetch

# Make models available
RUN mkdir -p /opt/mhcflurry2/models/
RUN cp -r /root/.local/share/mhcflurry/4/2.2.0/models_class1_presentation/models/* /opt/mhcflurry2/models/

# Clean up to reduce Docker image size
RUN apt-get clean
