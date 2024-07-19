FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    bzip2 \
    ca-certificates  \
    libpcre3 \
    libpcre3-dev \
    libssl-dev \
    zlib1g \
    zlib1g-dev \
    gfortran \
    openssl \
    make \
    build-essential \
    cmake \
    libncurses-dev \
    ncurses-dev \
    libbz2-dev \
    lzma-dev \
    liblzma-dev \
    libssl-dev \
    autoconf \
    automake \
    libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/walaj/svaba/releases/download/v1.2.0/svaba && \
    chmod +x svaba

ENV PATH /:$PATH
