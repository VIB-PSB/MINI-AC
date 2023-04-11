FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Nicolas MANOSALVA PEREZ"

# Install base utilities
4RUN apt-get update && \
5    apt-get install -y build-essentials  && \
6    apt-get install -y wget &&
7    apt-get clean && \
8    rm -rf /var/lib/apt/lists/*
9
10# Install miniconda
11ENV CONDA_DIR /opt/conda
12RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
13     /bin/bash ~/miniconda.sh -b -p /opt/conda
14
15# Put conda in path so we can use conda activate
16ENV PATH=$CONDA_DIR/bin:$PATH
