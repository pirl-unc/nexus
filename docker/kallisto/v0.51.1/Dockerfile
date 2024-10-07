# Use Rocky Linux 8 as the base image
FROM rockylinux:8

# Install dnf-plugins-core to get access to config-manager
RUN dnf install -y dnf-plugins-core

# Enable the powertools repository
RUN dnf config-manager --set-enabled powertools

# Install development tools and necessary libraries
RUN dnf groupinstall -y "Development Tools" && \
    dnf install -y \
        wget \
        bzip2 \
        ca-certificates \
        cmake \
        glibc-static \
        glibc-devel \
        zlib-devel \
        openssl-devel

# Download and Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Update PATH environment variable
ENV PATH=/opt/conda/bin:$PATH

# Update conda
RUN conda update -n base -c defaults conda

# Install kallisto dependencies
RUN conda config --set solver classic
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install -c conda-forge gcc_linux-64==14.1.0 hdf5==1.14.3
RUN conda install -c bioconda bustools==0.44.0

# Download, build, and install kallisto
RUN wget https://github.com/pachterlab/kallisto/archive/refs/tags/v0.51.1.tar.gz && \
    tar -xf v0.51.1.tar.gz && \
    cd kallisto-0.51.1 && \
    mkdir build && \
    cd build && \
    cmake .. -DENABLE_AVX2=OFF -DCOMPILATION_ARCH=OFF -DMAX_KMER_SIZE=64 && \
    make && \
    make install

# Clean up to reduce Docker image size
RUN conda clean --all
