FROM continuumio/miniconda:4.5.11
LABEL authors="Sheng-Da Hsu" \
      description="Docker image containing all requirements for freebayes-nf pipeline"
RUN apt-get update --fix-missing && apt-get install -y procps &&apt-get upgrade -y
COPY environment.yml /
ADD resources /tmp
RUN conda env create -f /environment.yml && \
    conda clean -a && \
    source activate gFB-1.0.2 && \
    cd /tmp/resources/fastq-stats-k && \
    make && \
    make install 
ENV PATH /opt/conda/envs/gFB-1.0.2/bin:$PATH
