# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

# Set noninteractive installation to avoid getting stuck during build
ARG DEBIAN_FRONTEND=noninteractive

# Install necessary packages for Miniconda installation
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends  \
    autoconf  \
    automake  \
    build-essential \
    bzip2 \
    ca-certificates \
    cmake  \
    curl \
    gfortran  \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev  \
    libncurses-dev \
    libpcre3 \
    libpcre3-dev \
    libssl-dev  \
    libssl-dev  \
    lzma-dev  \
    make  \
    ncurses-dev \
    openssl \
    wget \
    zlib1g  \
    zlib1g-dev  && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update conda

# Install R 3.6
RUN apt-get update && \
    apt-get install -y software-properties-common dirmngr && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" && \
    apt-get update && \
    apt-get install -y r-base=3.6* r-base-dev=3.6* && \
    apt-get clean

# Install Sequenza
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/cran/2020-04-01')); install.packages(c('readr', 'pbapply', 'squash', 'optparse', 'seqminer', 'BiocManager'), dependencies = TRUE)"

RUN Rscript -e "BiocManager::install('igordot/copynumber')"
RUN wget https://cran.r-project.org/src/contrib/Archive/sequenza/sequenza_3.0.0.tar.gz
RUN wget https://cran.r-project.org/src/contrib/Archive/iotools/iotools_0.3-1.tar.gz
RUN R CMD INSTALL iotools_0.3-1.tar.gz
RUN R CMD INSTALL sequenza_3.0.0.tar.gz
RUN conda install -c bioconda samtools==1.2
RUN conda install -c bioconda sequenza-utils==3.0.0

ADD run_sequenza.R /opt/sequenza/run_sequenza.R

# Clean up to reduce Docker image size
RUN apt-get clean
