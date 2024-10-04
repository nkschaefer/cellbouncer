FROM continuumio/miniconda3
RUN apt-get update && apt-get install -y build-essential zlib1g-dev libhts3 libhts-dev
COPY cellbouncer /cellbouncer
WORKDIR /cellbouncer
RUN conda env create --file=cellbouncer_minimum.yml
RUN make PREFIX=/usr/include
