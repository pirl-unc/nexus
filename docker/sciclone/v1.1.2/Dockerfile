FROM --platform=linux/amd64 ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# rgl (a sciClone dependency) is attached when library(sciClone) loads; in a
# headless container it must not try to open an X11/OpenGL display.
ENV RGL_USE_NULL=TRUE

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
      ca-certificates wget git \
      r-base r-base-dev \
      build-essential gfortran \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      zlib1g-dev libbz2-dev liblzma-dev \
      r-cran-ggplot2 r-cran-rgl r-cran-dplyr \
      r-cran-rcolorbrewer r-cran-plotrix r-cran-teachingdemos && \
    apt-get clean

RUN R -e 'install.packages(c("BiocManager", "remotes"), repos = "https://cloud.r-project.org")' && \
    R -e 'BiocManager::install(ask = FALSE, update = FALSE)' && \
    R -e 'BiocManager::install(c("IRanges", "limma"), ask = FALSE, update = FALSE)' && \
    R -e 'remotes::install_github("genome/bmm", upgrade = "never")'

RUN R -e 'options(repos = BiocManager::repositories()); remotes::install_url("https://github.com/genome/sciclone/archive/refs/tags/v1.1.2.tar.gz", dependencies = c("Depends", "Imports", "LinkingTo"), upgrade = "never")'

RUN R -e 'stopifnot(packageVersion("sciClone") == "1.1.2"); library(sciClone); cat("sciClone", as.character(packageVersion("sciClone")), "OK\n")'

RUN apt-get clean