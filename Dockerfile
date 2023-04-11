FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Nicolas MANOSALVA PEREZ"

# Install base utilities
RUN apt-get update
RUN apt-get install -y build-essentials wget
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH
