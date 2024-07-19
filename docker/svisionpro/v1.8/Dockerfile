# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 ca-certificates libgl1-mesa-glx libopengl0 \
    ffmpeg libsm6 libxext6 libpcre3 libgl1 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update -n base -c defaults conda

# Install SVision-pro
RUN pip install pysam==0.22.1
RUN pip install scipy==1.3.3
RUN pip install numpy==1.21.6
RUN pip install pillow==6.2.2
RUN pip install torch==1.10.2
RUN pip3 install opencv-python
RUN wget https://github.com/songbowang125/SVision-pro/archive/refs/tags/v1.8.tar.gz
RUN pip install v1.8.tar.gz
RUN tar -xzf v1.8.tar.gz

# Clean up to reduce Docker image size
RUN conda clean -all
RUN apt-get clean