# pantools-pipeline-v4
General purpose Snakemake pipeline for [PanTools](https://git.wur.nl/bioinformatics/pantools).

This pipeline can be used to run PanTools code for specific use cases and datasets.
This pipeline uses PanTools version 4.3.1. This pipeline is used for reproducibility of 
larger workflows in pangenomic analysis, and datasets provided by PanTools for tutorials and
reproducibility of experiments often come with their own configuration file.
More information is available below for setting up configuration for your own data and creating
a unique workflow.

Requirements: [Snakemake](https://snakemake.readthedocs.io/en/stable/), 
[Mamba](https://mamba.readthedocs.io/en/latest/).

## Installation
### Clone this git
For cloning this git, run:
```bash
git clone https://github.com/PanUtils/pantools-pipeline-v4
cd pantools-pipeline-v4
````

### Install Snakemake and Mamba
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

### all_variants
Create a pangenome with variation and run all major PanTools commands.

### panva
Create a pangenome and run all PanTools analysis functions required to make a 
[PanVA](https://github.com/PanBrowse/PanVA) instance.

## Creating your own workflow
It is advised not to use this pipeline the first time you create and analyse a pangenome with new data.
Often analysis during this process can alter the settings you want to use for different functions, and in-between 
analysis and parameter adjustments are not possible when using this automated workflow. The advised course of action
is to perform initial pangenome build and analysis steps without the pipeline, and if there is need for setting up a 
reproducible workflow for experimental replication or repeating the process with updated or similar data a custom 
configuration and workflow can be set up. You can store the commands options used in the initial experiment in your 
config, or otherwise find them stored in the log files of PanTools.
Important configuration parameters, and instructions for how to create your own workflow are discussed below. If you 
want to follow these steps, and keep up to date with the newest versions of PanTools, and this pipeline, we recommend 
setting up a [fork](
https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) of this 
repository. Inside of this fork, you can set up you own rules and standard configuration, while still being able to
pull in new developments if desired.

## Configuration
The configuration file contains four sections, which are further discussed below.

### Java settings
Settings for how to run the PanTools java program. The most important settings are the Xmx heap size and nice priority. 
Larger datasets (pangenomes with many and/or large genomes) will need a higher heap space to not run out of memory, 
especially for intensive processes such as the initial build step and homology grouping.
If your process runs out of memory, you need to increase the value for this parameter.
The [nice](https://ss64.com/bash/nice.html) setting lets you set a priority for the PanTools commands on a server, 
this is useful if you let large datasets run on a server with other traffic, so it can efficiently use more cores if 
there are more available without taking up all the server's processing power from other processes.
If you want to work with experimental, unreleased PanTools branches, this is possible by setting the jar file you want 
to run in the config, which will override the PanTools version in the conda environment used by the pipeline.

### Output
Three output directories can be set here. The results directory will contain your completed database. The construction 
directory can be the same as the results directory, in which case the pangenome will be constructed in one place.
If the construction directory differs from the results directory, the construction commands will be separated from 
analysis commands. Commands that write data to the Neo4j graph database are often transaction-heavy and work much better 
on a RAM or SSD disk, while for pangenome analysis functions there is no big difference. By setting separate 
directories, the construction commands will be performed on your RAM or SSD, after which the database will be moved 
over to the results directory for the analysis steps. This is recommended if space on your RAM or SSD is limited, or 
these disks are for temporary use only. 
Finally, a scratch directory can be set, which will be used for PanTools functions that create large amounts of 
temporary files. For larger datasets, the default temporary disk might not be large enough, so assigning a different 
directory is recommended.

### Parameters
In this section you can set the positional parameters available for PanTools, all of these are input files that are 
required for several PanTools functions. It is recommended to use full paths for your own files here.
These parameters are required for the following PanTools functions:
- **genomes**: [build_pangenome](https://pantools.readthedocs.io/en/stable/construction/build.html#build-pangenome) 
(required for all_pangenome, all_variants and panva)
- **annotations**: [add_annotations](https://pantools.readthedocs.io/en/stable/construction/annotate.html#add-annotations) 
(required for all_pangenome, all_variants and panva)
- **proteomes**: [build_panproteome](https://pantools.readthedocs.io/en/stable/construction/build.html#build-panproteome) 
(required for all_panproteome)
- **vcf**: [add_variants](https://pantools.readthedocs.io/en/stable/construction/annotate.html#add-variants) (optional)
- **pav**: [add_pav](https://pantools.readthedocs.io/en/stable/construction/annotate.html#add-pavs) (optional)
- **phenotypes**: [add_phenotypes](https://pantools.readthedocs.io/en/stable/construction/annotate.html#add-phenotypes) (optional)
- **phasing**: [add_phasing](https://pantools.readthedocs.io/en/stable/construction/annotate.html#add-phasing) (unused for preset workflows)
- **repeats**: [add_repeats](https://pantools.readthedocs.io/en/stable/construction/annotate.html#add-repeats) (unused for preset workflows )
- **synteny**: [add_synteny](https://pantools.readthedocs.io/en/stable/construction/synteny.html#add-synteny) (optional, unused for preset workflows)
- **blast**: [blast](https://pantools.readthedocs.io/en/stable/analysis/blast.html#blast) (unused for preset workflows)
- **short read (1 and 2)**: [map](https://pantools.readthedocs.io/en/stable/analysis/mapping.html#map) (required for all_pangenome and all_variants)
- **gene_selection**: Multiple PanTools commands have a ``--homology-file`` flag, which you can use to perform these 
  functions only on a selection of homology groups of interest. It is not currently possible to use gene identifiers
  instead of homology node IDs for homology selection of these commands, but this parameter can be used instead. 
  If you provide a file with comma-separated gene names for this parameter, the homology group identifiers for these 
  genes will be used for all functions that contain this option.

### Optional arguments
Here you can set all optional arguments (or flags) used for every function. With the exception of a few special 
arguments at the start of this section, they simply contain a string where you can put all command line options for any
function, and are named after every function (function_name.opts). For instance, if you want to run the 
gene_classification function with the following settings: \
``pantools gene_classification --unique-threshold=5 --core-threshold=95 tomato_DB``, the gene_classification.opts
parameter would be "--unique-threshold=5 --core-threshold=95".

The parameters that do not match this format are:
- **group.relaxation**: The value for the ``--relaxation`` option in pantools 
  [group](https://pantools.readthedocs.io/en/stable/construction/group.html#group). 
- **busco_protein.odb10**: The value for the ``--obd10`` option in pantools
  [busco_protein](https://pantools.readthedocs.io/en/stable/construction/group.html#busco-protein).
- **core_phylogeny.mode**: The value for the ``--clustering-mode`` option in pantools 
  [core_phylogeny](https://pantools.readthedocs.io/en/stable/analysis/phylogeny.html#core-phylogeny).

Either one of the group relaxation or BUSCO obd10 options is required. If the relaxation is not set, the pipeline will 
run busco_protein and optimal_grouping to determine the best relaxation setting for the provided database.
We recommend running busco_protein and optimal_grouping manually to determine a value for the grouping relaxation,
then setting the relaxation in the configuration file to make subsequent runs faster.
The core_phylogeny ``--clustering-mode`` option changes the output this command generates, which is why it must be set
separately from other core_phylogeny options. Make sure not to set any of these options in the list of optional 
parameters below.

## Custom workflows
You can create a custom workflow by adding a new rule to **workflow/Snakefile**. All rules that run a PanTools command 
create a ".done" file to signal their successful completion; a workflow can then be created simply by setting a list of 
done files for all PanTools functions you want to run in your workflow. The pipeline determines the order in which the
PanTools commands need to run, and commands that are required for other commands (like group for gene_classification), 
will automatically be added and do not need to be set. A .done file has the following naming convention:
``<results>/done/<database_type>.<command>.done``.
Here, **results** is the results directory set in the config file, **database_type** is either "pangenome" or 
"panproteome" and **command** is the name of the PanTools command you want to add to your workflow.
Use the existing rules as reference to create your own; it can then be run by adding the rule name to your Snakemake 
command.