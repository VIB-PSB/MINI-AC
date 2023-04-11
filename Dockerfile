FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Nicolas MANOSALVA PEREZ"

RUN apt-get update
RUN apt-get install -y wget make
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz
RUN tar -zxvf bedtools-2.27.0.tar.gz
WORKDIR bedtools2
RUN make
RUN make install
