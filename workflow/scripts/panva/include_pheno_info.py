#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Scripts:
# used to include phenotype/metadata columns for visualisation in PanVA.

# imports:
import pandas as pd
import os
from functools import reduce


# Functions:
def pheno_meta_maker(ptools_db_path):
    """
    Adds phenotype information of the (reference) genomes if available

    :param ptools_db_path: str - full path to the PanTools database
    :return: Dataframe - contains all meta and phenotype data of the genomes
    """
    pheno_over = os.path.join(ptools_db_path, 'phenotype_overview.txt')

    with open(pheno_over, 'r') as entry_line:
        pheno_lines = entry_line.readlines()
    entry_line.close()
    # find index for part 3
    # starts at one due to the blankline following the part3 header
    part3_index = 1
    for entry in pheno_lines:
        part3_index = part3_index + 1
        if entry.startswith('##Part 3'):
            break
    # subset only Part3 lines by found index
    part3_line_list = pheno_lines[part3_index:]

    # Extract the individual lines containing pheno info
    pheno_set_pattern = ": '"

    # For genome number and pheno node other string pattern
    pheno_num_pattern = "phenotype node:"

    # Add layer of security in the event Part3 in future won't be last part
    end_part_3_pattern = "##Part"

    # Dictionary
    pheno_dict = dict()
    pheno_dict['genome_nr'] = []
    pheno_dict['pheno_node_id'] = []

    for line in part3_line_list:
        # if new pheno id number is detected
        if pheno_num_pattern in line:
            genome_nr = line.split(',')[0].split(' ')[1]
            pheno_node = line.split(',')[1].strip('\n').split(' ')[-1]
            # add node_id and genome_nr to dictionary
            pheno_dict['genome_nr'].append(genome_nr)
            pheno_dict['pheno_node_id'].append(pheno_node)
        # if pheno feat id is detected
        elif pheno_set_pattern in line:
            # split the line in pheno_feat_id and the pheno feat code (MAYBE REMOVE spaces)
            pheno_feat_id = line.split(':')[0]
            pheno_feat_id = pheno_feat_id.lstrip()
            pheno_feat = line.split(':')[1].strip("\n '")
            # if the feature_id already exists add the feature code to the list
            if pheno_feat_id in pheno_dict:
                pheno_dict[pheno_feat_id].append(pheno_feat)
            # if the feature_id does not yet exist add it as a dict key and add
            else:
                pheno_dict[pheno_feat_id] = []
                pheno_dict[pheno_feat_id].append(pheno_feat)
        elif end_part_3_pattern in line:
            break
    # write the dictionary to a pandas dataframe
    df_phenos = pd.DataFrame.from_dict(pheno_dict)

    # typecast and json dump
    df_phenos['pheno_node_id'] = df_phenos['pheno_node_id'].astype(int)
    df_phenos = df_phenos.sort_values(by=['genome_nr'])

    return df_phenos


def acc_genome_metamerge(acc_meta_df, df_phenos):
    """
    Combines any (resequenced) accessions meta/ phenotype data with the other metadata/phenotype data.

    :param acc_meta_df: Dataframe - with at least one matching 'id' column and phenotype data
    :param df_phenos: Dataframe -  dataframe based on PanTools phenotype_overview standards.
    :return: Dataframe - combined version of the meta-/phenotype data
    """
    acc_meta_df_cols = acc_meta_df.columns
    pheno_cols = df_phenos.colums
    incommon_cols = [col for col in pheno_cols if col in acc_meta_df_cols]
    df_phenos = pd.merge(df_phenos, acc_meta_df, on=incommon_cols, how='outer')

    return df_phenos


# Ran in own pool
def hom_group_pheno(single_id_path, panva_path, df_phenos, df_seq_info):
    """
    Combines the individual sequence information with the corresponding meta/phenotype data.

    :param single_id_path: str - path to single homology grp id in alignments/../.. .
    :param panva_path: str - full path to output dir (serves as input for PanVA).
    :param df_phenos: Dataframe - containing all phenotype data
    :param df_seq_info: Dataframe - contains information on all the sequences
    :return: Dataframe & metadata.csv
    """
    # Get id
    hom_id = single_id_path.rsplit('/', 1)[1]
    # make output to panva
    hom_grp_pheno_pva = os.path.join(panva_path, hom_id)

    # merge information
    # hom_grp_pheno = hom_mrna_id_gnr.merge(df_phenos, on='genome_nr', how='left')
    hom_grp_pheno = pd.merge(df_seq_info[['genome_nr', 'mRNA_id']], df_phenos, on='genome_nr')

    hom_grp_pheno['genome_nr'] = hom_grp_pheno['genome_nr'].astype(int)

    hom_grp_pheno.to_csv(os.path.join(hom_grp_pheno_pva, 'metadata.csv'), index=False)

    return hom_grp_pheno
