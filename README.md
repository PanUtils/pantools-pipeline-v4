# pantools-pipeline-v4
General purpose Snakemake pipeline for PanTools v4.

This pipeline can be used to run PanTools code for specific use cases and datasets.

Requirements: Snakemake, Java, Mamba.

## Cloning this git
For cloning this git, run:
```bash
https://github.com/PanUtils/pantools-pipeline-v4
cd pantools-pipeline-v4
```

And check out a desired version (e.g. `v1.0.0`):
```bash
git checkout v1.0.0
```

## Install Snakemake and Mamba
If you don't have mamba, install it using
```bash
conda install -n base -c conda-forge mamba
```

Then, a Snakemake environment can be created using
```bash
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Which can be activated and verified with
```bash
conda activate snakemake
snakemake --help
```