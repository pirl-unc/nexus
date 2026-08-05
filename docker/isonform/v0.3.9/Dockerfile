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

RUN pip install isONform==0.3.9
RUN pip install isONclust==0.0.6.1
RUN pip install isONcorrect==0.1.3.5
RUN conda install -y --override-channels -c bioconda -c conda-forge spoa
ADD isON_pipeline.sh /isON_pipeline.sh
RUN chmod +x /isON_pipeline.sh

RUN conda install -y --override-channels -c conda-forge pandas==2.2.3
RUN pip install pysam==0.22.1

ADD run_isonform_clustered.py /run_isonform_clustered.py
RUN chmod +x /run_isonform_clustered.py

# Clean up to reduce Docker image size
RUN conda clean --all --yes
RUN apt-get clean