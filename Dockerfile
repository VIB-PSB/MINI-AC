FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Nicolas Manosalva Perez"

ADD requirements.txt requirements.txt

# Installing python 3.6 and pip3
RUN apt-get update

RUN apt-get install -y python3.6 python3-pip

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda install bedtools=2.27.0

#RUN apt-get install -y bedops=2.4.37 
#RUN apt-get install -y samtools=1.13 
#RUN pip3 install --upgrade pip

# Installing dependencies
RUN pip3 install -r requirements.txt
