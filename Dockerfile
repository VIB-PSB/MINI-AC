FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Nicolas Manosalva Perez"

ADD requirements.txt requirements.txt

# Installing python 3.6 and pip3
RUN apt-get update
# RUN apt-get update && \
#     apt-get install -y wget unzip gcc make zlib1g-dev && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*
RUN apt-get install -y python3.6 python3-pip
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz && \
    tar xvzf bedtools-2.27.0.tar.gz && \
    rm bedtools-2.27.0.tar.gz && \
    cd bedtools2 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bedtools2
#RUN apt-get install -y bedops=2.4.37 
#RUN apt-get install -y samtools=1.13 
#RUN pip3 install --upgrade pip

# Installing dependencies
RUN pip3 install -r requirements.txt
