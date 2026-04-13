FROM hkubal/clair3:v1.2.0

# Install Whatshap
RUN pip install --upgrade pip && \
    pip install pysam==0.23.0 && \
    pip install xopen>=1.2.0 && \
    pip install biopython==1.81 && \
    pip install whatshap==2.5

# Install samtools
RUN apt-get update && apt-get install -y libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xvjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18/ && \
    mkdir bin && \
    ./configure --prefix=$(pwd) && \
    make && \
    make install

# Add samtools to PATH
ENV PATH="/samtools-1.18:$PATH"

# Clean up to reduce Docker image size
RUN conda clean --all
RUN apt-get clean