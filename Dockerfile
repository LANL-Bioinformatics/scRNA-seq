# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

COPY environment.yml /

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc 


RUN conda env create -f environment.yml
RUN conda install -c conda-forge conda-pack


ADD scripts/Cluster_cell_type.py /opt/conda/envs/10x_covid_py/bin
ADD scripts/QC*.py /opt/conda/envs/10x_covid_py/bin

RUN conda-pack -n 10x_covid_py -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:buster AS runtime

COPY --from=build /venv /venv
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash