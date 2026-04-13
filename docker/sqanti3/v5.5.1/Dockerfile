FROM anaconesalab/sqanti3:v5.5.1

USER root

RUN conda install -n base conda-libmamba-solver \
    && conda config --set solver libmamba

COPY SQANTI3.conda_env.yml /tmp/SQANTI3.conda_env.yml

RUN conda env update -n sqanti3 -f /tmp/SQANTI3.conda_env.yml \
    && conda clean -afy

# Activate sqanti3 env by default so singularity exec works
ENV PATH=/conda/miniconda3/envs/sqanti3/bin:$PATH
ENV CONDA_DEFAULT_ENV=sqanti3
ENV CONDA_PREFIX=/conda/miniconda3/envs/sqanti3

ENTRYPOINT []

USER user
