# Use Ubuntu 18.04 LTS as the base image
FROM ubuntu:18.04

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

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# Test conda installation by updating conda
RUN conda update conda

# Install R
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.3.tar.gz && \
    tar -xf R-3.6.3.tar.gz && \
    cd R-3.6.3/ && \
    ./configure \
        --prefix=/opt/R/R-3.6.3 \
        --with-readline=no \
        --with-x=no && \
    make && \
    make install

# Add R to PATH
ENV PATH /opt/R/R-3.6.3/bin:$PATH

# Install sequenza
RUN Rscript -e "install.packages('pbapply', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('squash', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('vroom', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('readr', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('optparse', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('seqminer', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
RUN Rscript -e "install.packages('BiocManager', repos = 'http://cloud.r-project.org', dependencies = TRUE)"
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
