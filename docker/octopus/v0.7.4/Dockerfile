FROM dancooke/octopus:invitae--eae1ab48_0

ENV DEBIAN_FRONTEND=noninteractive

# Point EOL Ubuntu repos to old-releases
RUN sed -i \
      -e 's|http://archive.ubuntu.com/ubuntu|http://old-releases.ubuntu.com/ubuntu|g' \
      -e 's|http://security.ubuntu.com/ubuntu|http://old-releases.ubuntu.com/ubuntu|g' \
      /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN samtools --version

ENTRYPOINT []
