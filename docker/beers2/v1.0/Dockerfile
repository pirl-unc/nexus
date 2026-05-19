FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 unzip ca-certificates  \
    libpcre3 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_26.3.2-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

# BEERS2 v1.0 hard-pins pysam==0.20.0, whose cy_build.py does
# `from pkg_resources import Distribution`. setuptools >= 81 (late 2025)
# dropped the vendored pkg_resources from clean installs, so pip's
# build-isolation env (which pulls the latest setuptools by default)
# fails to build pysam. Constrain setuptools in the build env so
# pkg_resources stays importable. PIP_CONSTRAINT is honored by both
# the outer install and the isolated PEP 517 build envs pip spins up
# per-package.
RUN printf 'setuptools<81\nwheel\n' > /etc/pip-constraints.txt
ENV PIP_CONSTRAINT=/etc/pip-constraints.txt

# Download and unpack BEERS2
RUN wget --quiet https://github.com/itmat/BEERS2/archive/refs/tags/v1.0.tar.gz && \
    tar -xvf v1.0.tar.gz && \
    cd BEERS2-1.0/ && \
    pip install .

RUN conda clean -a -y && \
    apt-get clean
