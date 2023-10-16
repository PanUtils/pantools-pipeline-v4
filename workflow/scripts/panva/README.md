# PanToVA 2022 2023
##### Sander Vlugter

WUR - Bioinformatics master thesis 2022-2023.

### Functionality
Collection of scripts used to pre-process (PanTools) pangenome data into PanVA input data.

### How to use
Construct/activate the conda environment provided in envs/ of PanToVA. Copy the template config file from the 
config_template folder to a directory of choice and fill in the config based on your pangenome data. Examples of the
config file are also provided within the config_template folder.
If any setting in the config file are unclear please consult the "config help" document under docs (not yet available).

After activating the conda environment and filling out the config you can start the pre-processing by using:

``
python3 path/to/pan_to_va.py path/to/config.ini
``

### Disclaimers:

Some of the tests sets have older file versions which can cause errors. If users want to use these older
versions, you need to manually indicate the version of the files in the scripts. 
Support for the older PanTools file versions is going to be removed in future..
