#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:

# Imports:
from homgroup_seq_info import *
import pheno_specific_var
from make_alignments import *
from include_pheno_info import *


# ran in pool
def prep_group_pheno(hom_id, panva_p, seqtype='nuc_trimmed', pheno_var=None, meta=None):
    """
    Flashes out the alignments for a homology group, by adding additional information if present.

    :param hom_id: str - path to single homology grp id in alignments/../.. .
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :param seqtype: str - define the type of sequence default 'nuc_trimmed'.
    :param pheno_var: bool - indicates if there is data on meta/phenotype specific variant positions.
    :param meta: optional - can be None or dataframe containing meta/phenotype data on the aligned sequences.
    :return: meta_info - Dataframe - basis of the homologies.json
    """

    df_seq_info = create_seq_info(hom_id, panva_p)

    # makes sequences.csv
    df_all_seq = merge_sequences(hom_id, panva_p, method=seqtype)

    all_info_seq = pd.merge(df_all_seq, df_seq_info, on=['mRNA_id'])

    # if both not none ADD phenometa data AND pheno specific var
    if meta is not None and pheno_var is not None:
        # makes metadata.csv
        # print("prep_group: meta present and pheno_var present")
        new_thing = hom_group_pheno(hom_id, panva_p, meta, df_seq_info)
        meta_info = alignment_posinfo(hom_id, panva_p, all_info_seq, method=seqtype, phe_var=pheno_var, meta=new_thing)

    # if meta not None add phenometa but not pheno var
    elif meta is not None and pheno_var is None:
        # print("prep_group: meta present and pheno_var NOT present")
        # should work otherwise make if statement based on pheno_var
        new_thing = hom_group_pheno(hom_id, panva_p, meta, df_seq_info)
        meta_info = alignment_posinfo(hom_id, panva_p, all_info_seq, method=seqtype, phe_var=pheno_var, meta=new_thing)

    # if meta is None and pheno var is None just do positions
    else:
        # print("prep_group: meta NOT present and pheno_var NOT present")
        meta_info = alignment_posinfo(hom_id, panva_p, all_info_seq, method=seqtype, phe_var=pheno_var, meta=meta)

    return meta_info
