#!/usr/bin/python

# Sander Vlugter
# Date start: 21-12-2022
# WUR Bioinformatics CPG

# Imports
import os
from pheno_specific_var import *
from include_var_pos import *
import pandas as pd


# Functions
def alignment_posinfo(single_id, panva_p, all_seq_all_info, method='nuc_trimmed', phe_var=None, meta=None):
    """
    Adds phenotype/metadata information to the positions in the alignment if available and parse alignment data.

    :param single_id: str - path to single homology grp id in alignments/../.. .
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :param all_seq_all_info: Dataframe - combination of sequence and if exists meta information
    :param method: str - define the type of sequence default 'nuc_trimmed'.
    :param phe_var: bool - indicates if there is data on meta/phenotype specific variant positions.
    :param meta: optional - can be None or dataframe containing meta/phenotype data on the aligned sequences.

    :return:
    """
    # keep track of hom_id
    # hom_id = single_id.rsplit('/', 1)[1]
    # print(hom_id)
    # output dir
    # panva_out = os.path.join(panva_p, hom_id)

    # list to hold temp dataframes
    holding_tmp_df = []

    method = method + "_seq"

    # all_seq_all_info = pd.merge(df_all_seq, df_seq_info, on=['mRNA_id'])
    # print(all_seq_all_info)
    for mrna_id in all_seq_all_info['mRNA_id']:
        # get genome number
        gnme_nr = all_seq_all_info[all_seq_all_info['mRNA_id'] == mrna_id]['genome_nr'].values[0]
        # print(gnme_nr)
        # get sequence as nuc per elem list
        if method == 'var_trimmed_seq':
            method = 'nuc_trimmed_seq'
        sequence = [*all_seq_all_info[all_seq_all_info['mRNA_id'] == mrna_id][method].to_list()[0]]
        # print(sequence)
        # make df
        all_pos_df_tmp = pd.DataFrame(sequence, columns=['nucleotide'])
        # print(all_pos_df_tmp)
        # grab index and add 1 (position in alignment)

        all_pos_df_tmp['position'] = all_pos_df_tmp.index + 1
        # reset index
        # all_pos_df_tmp = all_pos_df_tmp.reset_index()
        # print(all_pos_df_tmp)
        # rename the 'index' column to position
        # rename the 'index'column to position
        # all_pos_df_tmp.rename(columns={"index": 'position'})
        all_pos_df_tmp['mRNA_id'] = mrna_id
        all_pos_df_tmp['genome_nr'] = gnme_nr
        # print(all_pos_df_tmp)
        holding_tmp_df.append(all_pos_df_tmp)

    # combine them
    all_pos_df_unorder = pd.concat(holding_tmp_df)
    # change column order to match
    all_pos_df = all_pos_df_unorder.reindex(columns=['mRNA_id', 'genome_nr', 'position', 'nucleotide'])
    # maybe:
    # all_pos_df['genome_nr'] = all_pos_df['genome_nr'].astype(int)
    # print('after make alignments')
    # print(all_pos_df)
    # print(phe_var)
    if meta is not None and phe_var is not None:
        # print("meta present phe not none")
        all_pos_df = add_pheno_specificity(single_id, all_pos_df, meta, pheno_var=phe_var)
        # found in include_var_pos.py
        meta_info = create_var_pos_df(single_id, panva_p, all_pos_df, pheno=phe_var, version='new')
        # all_pos_df.to_csv(os.path.join(panva_out, 'alignments.csv'), index=False)
    elif meta is not None and phe_var is None:
        # print("meta present phe is none")
        meta_info = create_var_pos_df(single_id, panva_p, all_pos_df, pheno=phe_var, version='old')
        # all_pos_df.to_csv(os.path.join(panva_out, 'alignments.csv'), index=False)
    elif meta is None and phe_var is None:
        # print("both none")
        meta_info = create_var_pos_df(single_id, panva_p, all_pos_df, pheno=phe_var, version='old')
        # all_pos_df.to_csv(os.path.join(panva_out, 'alignments.csv'), index=False)
    return meta_info
