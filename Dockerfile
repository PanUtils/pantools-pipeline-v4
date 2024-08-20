FROM condaforge/mambaforge:23.11.0-0

RUN apt-get update -y
RUN apt-get install -y git
RUN conda install -c conda-forge -y openjdk

RUN git clone https://github.com/PanUtils/pantools-pipeline-v4.git
WORKDIR pantools-pipeline-v4

RUN mamba create -c conda-forge -c bioconda -n snakemake pantools snakemake
RUN echo "conda activate snakemake" >> ~/.bashrc

