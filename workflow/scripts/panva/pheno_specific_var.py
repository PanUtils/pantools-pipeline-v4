#!/usr/bin/python

# Sander Vlugter
# Date start: 21-12-2022
# WUR Bioinformatics CPG

# Script in preparation of the updated phenotype specific snp formats. (subject to change)


# Imports:

import os
import pandas as pd
import numpy as np
import logging
from functools import reduce
# Function(s):


def add_pheno_specificity(single_id_path, al_pos, phenos_df, pheno_var=None):
    """

    :param single_id_path:
    :param al_pos:
    :param phenos_df:
    :param pheno_var:
    :return:
    """
    # path to pantools_db/id/output/ location of the {pheno}_spec_var
    hom_id_output = os.path.join(single_id_path, "output/phenotype/")

    # place to filter specific phenotypes
    if isinstance(pheno_var, str):
        if pheno_var == 'all' or pheno_var == ' all':
            # create list of all selected groups
            pheno_var_file_id = [file for file in os.listdir(hom_id_output) if os.path.isfile(os.path.join(hom_id_output, file))]
        elif os.path.isfile(pheno_var):
            with open(selection, 'r') as grp_file:
                pheno_var_file_id = [line.strip()for line in grp_file]
            grp_file.close()
            # check if files were given otherwise try to deal with it.
        else:
            logging.error("The pheno var in the config file can be empty 'all' or '.txt listing specific pheno files'.")
            raise ValueError("The pheno var in the config empty, 'all' or '.txt listing the specific pheno files'.")

    # place a check if the phenotype from file name matches the phenotypes in the pheno_df
    phenovar_fullpath = []
    for pheno in pheno_var_file_id:
        #print(pheno)
        #pheno_out_folder = os.path.join(hom_id_output)
        path_to_pheno_out = os.path.join(hom_id_output, pheno)
        #print(path_to_pheno_out)
        phenovar_fullpath.append(os.path.join(path_to_pheno_out))

    # merge aligned pos information with phenotype information

    align_pos_pheno = pd.merge(al_pos, phenos_df, on=['mRNA_id', 'genome_nr'], how='outer')

    # print("aligned with pheno mate \n", align_pos_pheno)

    # list holding all phenotype variant pos classifiers
    pheno_holding_list = []

    # change it to the list phenovar_fullpath
    seq_type = 'nuc_trimmed'
    only_nuc_files = []
    for file in phenovar_fullpath:
        only_nuc_files.append(file)
    #for file in phenovar_fullpath:
    #    if seq_type in file:
    #        only_nuc_files.append(file)

    # for all files check which match the file naming scheme.
    for file in only_nuc_files:
        pheno_id = file.rsplit("/", 1)[1]
        # keeping only the phenotype key word
        if file.startswith("nuc_"):
            pheno_id = pheno_id.split("nuc_trimmed_")[1]
        elif file.startswith("var_"):
            pheno_id = pheno_id.split("var_trimmed_")[1]
        pheno_id = pheno_id.rsplit("_", 1)[0]

        # make column indicator from the pheno_id
        col_pheno = pheno_id + "_var"
        # full path and open the .csv
        specific_pheno_var = pd.read_csv(file)
        specific_pheno_var.rename(columns={'Position': 'position', 'Letter': 'nucleotide', 'Specific': 'specific', 'Exclusive': 'exclusive'}, inplace=True)
        # select only the required columns & merge with the alignment
        spec_pheno_align = pd.merge(align_pos_pheno, specific_pheno_var[['position', 'nucleotide', 'specific', 'exclusive']],
                                    on=['position', 'nucleotide'])
        # change the elements for col 'specific' as definition makes other values impossible
        spec_pheno_align.loc[spec_pheno_align['specific'].notnull(), 'specific'] = 'specific'
        # evaluate if the phenotype category (of the aligned sequence position) == the cat in exclusive
        spec_pheno_align.eval('exclusive = @pheno_id == exclusive', inplace=True)
        # replace Bool with tag 'exclusive' for matching pheno_cat pairs
        spec_pheno_align['exclusive'] = spec_pheno_align['exclusive'].map({True: 'exclusive', False: None})
        # add new pheno_col indicating category of the nucleotide in relation to phenotype
        spec_pheno_align[col_pheno] = spec_pheno_align['specific'].fillna(spec_pheno_align['exclusive'])

        # if required option to change None fill here
        # print("yes this is spec_align")
        # print(spec_pheno_align)
        # drop the temp columns
        spec_pheno_align = spec_pheno_align.drop(columns=['exclusive', 'specific', pheno_id], axis=1, errors='ignore')

        # deposit the df in the holding list for storing till all phenotypes have been assessed.
        pheno_holding_list.append(spec_pheno_align)


    accepted_columns = ['mRNA_id', 'genome_nr', 'nucleotide', 'position']
    df_phenos_columns = phenos_df.columns
    remove_cols = [col for col in df_phenos_columns if col not in accepted_columns]
    pheno_hold_2 = []
    for hom_grp in pheno_holding_list:
        hom_grp = hom_grp.drop(columns=remove_cols, errors='ignore')
        pheno_hold_2.append(hom_grp)


    # combine the individual phenotype dataframes
    al_pos_all_pheno = reduce(lambda df1, df2: pd.merge(df1, df2, how='outer'), pheno_hold_2)
    all_pos_all_pheno = pd.merge(al_pos, al_pos_all_pheno, how='outer')
    # print("actual all pos with pheno aligned /n", all_pos_all_pheno)
    # return to make alignments
    return all_pos_all_pheno
