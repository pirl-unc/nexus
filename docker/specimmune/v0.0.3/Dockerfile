FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends wget curl git bzip2 unzip ca-certificates  \
    libpcre3 libpcre3-dev libssl-dev zlib1g zlib1g-dev gfortran openssl \
    make build-essential cmake libncurses-dev ncurses-dev \
    libbz2-dev lzma-dev liblzma-dev libssl-dev autoconf automake libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_24.1.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/lib:$PATH

RUN conda update -n base -c defaults conda

# Download and unpack SpecImmune v0.0.3 release
WORKDIR /opt
RUN wget --quiet https://github.com/deepomicslab/SpecImmune/archive/refs/tags/v0.0.3.tar.gz && \
    tar -xvf v0.0.3.tar.gz && \
    rm v0.0.3.tar.gz

# Create the SpecImmune conda environment from the bundled environment.yml
# (Mirrors install.sh: `conda env create -n SpecImmune -f environment.yml`.)
RUN conda env create -n SpecImmune -f /opt/SpecImmune-0.0.3/environment.yml

# Install dysgu into the SpecImmune env without pulling in its own deps
# (deps are already resolved by environment.yml).
RUN conda run -n SpecImmune pip install --no-deps dysgu==1.6.2

# Activate the SpecImmune env by default for any downstream RUN / docker exec
ENV PATH=/opt/conda/envs/SpecImmune/bin:$PATH
ENV CONDA_DEFAULT_ENV=SpecImmune
ENV CONDA_PREFIX=/opt/conda/envs/SpecImmune

WORKDIR /opt/SpecImmune-0.0.3

# Build the HLA database into /opt/SpecImmune-0.0.3/db/ so the image ships
# self-contained (no first-run setup required).
RUN python scripts/make_db.py -o ./db -i HLA
RUN python scripts/make_db.py -o ./db -i KIR
RUN python scripts/make_db.py -o ./db -i IG_TR

RUN conda clean -a -y && \
    apt-get clean