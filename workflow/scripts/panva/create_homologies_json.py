#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:
# Building the homologies json required in PanVA

# imports:
import pandas as pd
import numpy as np
import os
import json
from functools import reduce
# Functions:


def gene_classification(gene_class_path, filtered_id_list):
    """
    Collect pangenome classification labels for each of the homology groups.

    :param gene_class_path: full path to the method of gene classification of choice.
    :param filtered_id_list: list of the homology group ids
    :return: Dataframe - containing id and the 'class'ification of the homology group
    """
    if os.path.isfile(gene_class_path):
        classified_df = pd.read_csv(gene_class_path)

    else:
        raise ValueError("The problem in path to the 'gene classification'")
    clean_id_list = []
    for hom_id in filtered_id_list:
        hom_id = hom_id.rsplit('/', 1)[1]
        hom_id = hom_id.strip()
        clean_id_list.append(hom_id)
    clean_id_list = list(map(int, clean_id_list))

    classified_df.columns = classified_df.columns.str.replace(' ', '_').str.lower()
    # print(classified_df.columns)
    hom_classified = classified_df[classified_df['homology_group_id'].isin(clean_id_list)][['homology_group_id', 'class']]

    hom_classified.columns = ['id', 'classification']
    hom_classified['id'] = hom_classified['id'].astype(str)

    # print(hom_classified)
    return hom_classified


def create_homologies(pangenome_path, panva_p, df_aligned_filt, hom_classified, pheno_var=None):
    """
    Building the homologies.json required in PanVA

    :param pangenome_path: str - path to the pangenome db
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :param df_aligned_filt: list - containing the homology group ids that passed filtering
    :param hom_classified: Dataframe - dataframe produced by gene_classification
    :param pheno_var: None or List - if list expects aligned pos information on specific or exclusive
    :return: homologies.json
    """
    group_inf_path = os.path.join(pangenome_path, "group_info/group_info.txt")

    hom_ids = df_aligned_filt.id.values.tolist()

    with open(group_inf_path, 'r') as grp_inf:
        info_lines = grp_inf.readlines()
    grp_inf.close()

    # hold the different information of interest
    grp_id_list = []
    gene_name_list = []
    function_list = []
    # num_mem_list = []
    # num_genomes_list = []
    # protein length could still be added

    # extraction of information indicator
    extract_bool = False

    for line in info_lines:
        if line.startswith("#Homology"):
            # check if the homology group is in the selected groups
            checker = line.rsplit(" ", 1)[1]
            checker = checker.strip()
            if checker in hom_ids:
                extract_bool = True
        if extract_bool:
            if line.startswith("#Homology"):
                grp_id = line.rsplit(" ", 1)[1]
                grp_id = grp_id.strip()
                grp_id_list.append(grp_id)
            if line.startswith("All gene names:"):
                gene_name = line.split("names:")[1]
                gene_name = gene_name.strip('\n')
                gene_name_list.append(gene_name)
            if line.startswith("All functions:"):
                all_func = line.split(": ", 1)[1]
                all_func = all_func.strip('\n')
                if all_func == '':
                    all_func = '-'
                function_list.append(all_func)
            # if line.startswith('num_members:'):
            #     # includes multi copy appearances aka should be called something like: number aligned seq
            #     num_mem = int(line.split(": ", 1)[1])
            #     num_mem_list.append(num_mem)
            # if line.startswith("Found in"):
            #     # get the back of the line
            #     num_genome_members = line.split("in ", 1)[1]
            #     # get just the number
            #     num_genome_members = int(num_genome_members.split("genomes", 1)[0])
            #     # print(num_genome_members)
            #     num_genomes_list.append(num_genome_members)
        if line.startswith('\n'):
            extract_bool = False
    # data_json = {"id": grp_id_list, "name": gene_name_list, "members": num_mem_list,
    #              "all_functions": function_list, "in_num_genomes": num_genomes_list}
    data_json = {"id": grp_id_list, "name": gene_name_list, "all_functions": function_list}

    homologies_df = pd.DataFrame(data_json)
    df_aligned_filt.rename(columns={'homology_id': 'id'})

    bighomologies_df = pd.merge(homologies_df, df_aligned_filt, on='id')

    base_info = pd.merge(bighomologies_df, hom_classified, on='id')
    base_info = base_info.rename(columns={'numb_members': 'members', 'align_length': 'alignment_length'})

    # print(pheno_var)
    if pheno_var is not None:
        pheno_var_df = reduce(lambda df_1, df_2: pd.merge(df_1, df_2, how='outer'), pheno_var)
        # pheno_var_df = pheno_var_df.fillna(False)
        # print(pheno_var_df)
        base_info = pd.merge(base_info, pheno_var_df, on='id')

    first_col = ['id', 'alignment_length', 'members']
    other_cols = [col for col in base_info.columns if col not in first_col]
    base_info = base_info[first_col + other_cols]
    # print(base_info)

    hom_json = base_info.groupby(['id', 'alignment_length', 'members']).apply(
        lambda x: x.iloc[:, 3:].to_dict('records')[0]).reset_index().rename(columns={0: 'metadata'}).to_json(
        orient='records')

    base_info.to_csv(os.path.join(panva_p, 'base_info.csv'), index=False)

    with open(os.path.join(panva_p, 'homologies.json'), 'w', encoding='utf-8') as new_json:
        json.dump(json.loads(hom_json), new_json, indent=2, ensure_ascii=False)
    # homologies_df = pd.merge(homologies_df, homid_align_len, on='id')
    # homologies_json = homologies_df.to_dict('records')
    # print(homologies_json)
    return base_info
