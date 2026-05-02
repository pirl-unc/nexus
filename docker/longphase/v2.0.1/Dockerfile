# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Minimal runtime + tools needed to fetch & extract the longphase release.
# tabix provides bgzip + tabix (used to compress/index the phased VCF).
# samtools is used to index the haplotagged BAM and extract a (read, HP) TSV.
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        wget \
        ca-certificates \
        xz-utils \
        tabix \
        samtools \
        gawk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install longphase (prebuilt x86_64 binary)
RUN wget --quiet https://github.com/twolinin/longphase/releases/download/v2.0.1/longphase_linux-x64.tar.xz \
        -O /tmp/longphase.tar.xz && \
    tar -xJf /tmp/longphase.tar.xz -C /usr/local/bin/ && \
    chmod +x /usr/local/bin/longphase_linux-x64 && \
    ln -s /usr/local/bin/longphase_linux-x64 /usr/local/bin/longphase

# Sanity check
RUN longphase --help > /dev/null 2>&1 || longphase 2>&1 | head -1
