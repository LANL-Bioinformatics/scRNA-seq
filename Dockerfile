# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

COPY environment_v3.yml /

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc 


RUN conda env create -f environment_v3.yml

# #2 packages have no conda versions, installing into environment with pip
# RUN conda install pip
# RUN /opt/conda/envs/10x_covid_py/bin/pip install pydeseq2==0.4.11
# RUN /opt/conda/envs/10x_covid_py/bin/pip install pertpy==0.9.4
# #fix pip clobbering files
# RUN conda install -n 10x_covid_py -c conda-forge numpy=1.26.4 --force-reinstall
# RUN conda install -n 10x_covid_py -c conda-forge numpy-base=1.26.4 --force-reinstall
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
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash