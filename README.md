# pantools-pipeline-v4
General purpose Snakemake pipeline for PanTools v4.

This pipeline can be used to run PanTools code for specific use cases and datasets.

Requirements: Snakemake, Mamba.

## Cloning this git
For cloning this git, run:
```bash
git clone https://github.com/PanUtils/pantools-pipeline-v4
cd pantools-pipeline-v4
````

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

## Running the pipeline
## Run the pipeline
The pipeline can be run with

```bash
snakemake [rule] --use-conda --cores <threads> [--configfile <config>] [--conda-prefix <prefix>] [--until <function>]
```

**threads**: number of threads to use. \
**config**: custom configuration file. \
**prefix**: path to directory containing conda environments for reuse. \
**function**: name of a pantools function, the pipeline stops if this function is complete. \

If no config is provided, the pipeline will run on a small yeast test dataset.
The possible rules are discussed below. The pipeline will run all major PanTools functions 
for a pangenome if no rules are specified.

## Rules
### all_pangenome
Create a pangenome and run all major PanTools analysis functions for pangenomes.

### all_panproteome
Create a panproteome and run all major PanTools analysis functions for panproteomes.

### panva
Create a pangenome and run all PanTools analysis functions required to make a PanVa instance.