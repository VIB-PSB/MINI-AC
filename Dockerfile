FROM docker.io/mambaorg/micromamba@sha256:cc777a148e0de55b9a44677464d195c01a77e1a9ef0f047416b57f7b6d95029c

ADD conda_mini_ac.yml conda_mini_ac.yml

USER root

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps wget
RUN apt-get install -y python3.6 

# Installing dependencies
RUN micromamba create -f conda_mini_ac.yml && \
micromamba clean -a

ENV PATH /opt/conda/envs/mini_ac/bin:$PATH

