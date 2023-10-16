#!/usr/bin/python

# Sander Vlugter
# Date start: 21-12-2022
# WUR Bioinformatics CPG
# Function(s):
# creates the variables.csv file combining data of the alignment file
# with the variant inf positions file.

# Remarks:
# Optionally can add the variable position information as col to alignment.

# Imports:
import os
import pandas as pd
import shutil
from collections import defaultdict
# Functions:


def create_var_pos_df(single_id, panva_p, df_all_pos, pheno=None, version='old'):

    # identify single id
    hom_id = single_id.rsplit('/', 1)[1]
    # path to output
    hom_id_out = os.path.join(panva_p, hom_id)
    # path to data
    # had to do it this way because some versions used old naming convention (nuc vs nucleotide)
    var_df_path = os.path.join(single_id, "output/var_inf_positions/nuc_trimmed_variable_positions.csv")
    if not os.path.isfile(var_df_path):
        var_df_path = os.path.join(single_id, "output/var_inf_positions/nucleotide_trimmed_variable_positions.csv")
    if not os.path.isfile(var_df_path):
        var_df_path = os.path.join(single_id, "output/var_inf_positions/variants_trimmed_variable_positions.csv")
    meta_info = []
    #print(single_id)
    if version == 'old' and pheno is None:
        #df_var_pos = pd.read_csv(var_df_path)
        if not os.path.isfile(var_df_path):
            print("the file nuc_trimmed_variable_pos... does not exist for{}".format(hom_id))
        df_var_pos = pd.read_csv(var_df_path)
        # print(hom_id)

        if "Informative" not in df_var_pos.columns:
            #print(hom_id)
            df_all_pos.to_csv(os.path.join(hom_id_out, 'alignments.csv'), index=False)
        else:
            df_var_pos = df_var_pos.rename(columns={"Position": "position", "Informative": "informative"})
            df_var_pos['informative'] = [True if i == 'I' else False for i in df_var_pos['informative']]
            df_var_pos.to_csv(os.path.join(hom_id_out, 'variable.csv'), index=False)
            df_all_pos.to_csv(os.path.join(hom_id_out, 'alignments.csv'), index=False)
        meta_info = [hom_id]

    elif version == 'new' and pheno is None:
        new_name = os.path.join(hom_id_out, "variable.csv")
        shutil.copyfile(var_df_path, new_name)
        df_all_pos.to_csv(os.path.join(hom_id_out, 'alignments.csv'), index=False)
        meta_info = [hom_id]

    elif version == 'new' and pheno is not None:
        df_var_pos = pd.read_csv(var_df_path)
        df_var_pos = df_var_pos.rename(columns={"Position": "position", "Informative": "informative"})
        var_cols = df_all_pos.columns
        static_cols = ['mRNA_id', 'genome_nr', 'position', 'nucleotide']
        pheno_cols = [col for col in var_cols if col not in static_cols]
        meta_info = [("id", hom_id)]
        #print(df_var_pos)
        # print("mate all pos below \n", df_all_pos)
        #print(df_all_pos)
        for pheno in pheno_cols:
            # hard to read but relatively fast
            df_var_pos[pheno] = df_var_pos['position'].isin(df_all_pos[~df_all_pos[pheno].isnull()]['position'])
            if True in df_var_pos[pheno].unique():
                meta_info.append((pheno, True))
            else:
                meta_info.append((pheno, False))
        #meta_info = pd.DataFrame({pheno: val} for pheno, val in meta_info).bfill().reset_index(drop=True)
        dict_it = defaultdict(list)
        for col, val in meta_info:
            dict_it[col].append(val)
        meta_info = pd.DataFrame(dict_it)

        df_var_pos.to_csv(os.path.join(hom_id_out, 'variable.csv'), index=False)
        df_all_pos.to_csv(os.path.join(hom_id_out, 'alignments.csv'), index=False)
        # print('this is meta_info')
        # print(meta_info)
    return meta_info
