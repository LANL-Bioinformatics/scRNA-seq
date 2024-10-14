# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

COPY environment_v3.yml /

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc 


RUN conda env create -f environment_v3.yml

#create packed environment with scripts
RUN conda install -c conda-forge conda-pack


ADD scripts/* /opt/conda/envs/10x_sc/bin

RUN conda-pack -n 10x_sc -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

#runtime stage to compress layers+improve image size
FROM debian:latest AS runtime

COPY --from=build /venv /venv
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash