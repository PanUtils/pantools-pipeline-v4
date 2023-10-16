#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Scripts:
# collect and report alignment information from the pre-filtered selection of groups.

# Imports:
import os
import pandas as pd
import logging

# Functions:


def get_alignment_info(single_id_path, seqtype='nuc_trimmed'):
    """
    Gathers alignment info for building the homologies.json and to filter on alignment statistics later.
    *Used in pool to speed up the process.

    :param single_id_path: full path to a single homology group in the alignments extension
    :return: list of the alignment stats for the homology group
    """
    # keep just the id part
    single_id = single_id_path.rsplit('/', 1)[1]
    # add the correct extension to output/nuc_trimmed_alignment.info
    if seqtype == 'nuc_trimmed':
        single_id_align_info = os.path.join(single_id_path, "output/nuc_trimmed_alignment.info")
        if not os.path.exists(single_id_align_info):
            single_id_align_info = os.path.join(single_id_path, "output/nucleotide_trimmed_alignment.info")
    elif seqtype == 'var_trimmed':
        single_id_align_info = os.path.join(single_id_path, "output/variants_trimmed_alignment.info")
    else:
        raise ValueError("Issue in finding the correct '[]_alignment.info' file. nuc_trimmed or variants_trimmed")

    with open(single_id_align_info, 'r') as align_info:
        informative_lines = [next(align_info) for x in range(0, 5)]
    align_info.close()
    # extract information for homologies.json
    num_seq = informative_lines[0].strip()
    num_seq = int(num_seq.split(':')[1])

    num_genomes = informative_lines[1].strip()
    num_genomes = int(num_genomes.split(':')[1])

    align_len = informative_lines[2].strip()
    align_len = int(align_len.split(':')[1])

    num_var_pos = informative_lines[3].strip()
    num_var_pos = int(num_var_pos.split(':')[1])

    num_inf_pos = informative_lines[4].strip()
    num_inf_pos = int(num_inf_pos.split(':')[1])

    combined_align_info = [single_id, num_seq, num_genomes, align_len, num_var_pos, num_inf_pos]
    return combined_align_info


def filter_on_alignment(combined_align_info, pan_grp_path, panva_path, align_len=90, num_members=2, num_uni_mem=1):
    """
    Filters the selected homology groups on different alignment based criteria

    :param combined_align_info: list - of lists of alignment stats
    :param pan_grp_path: str - full path to alignment/<msa_type>/grouping
    :param panva_path: str - full path to the output dir (serves as input for PanVA)
    :param align_len: int - value minimal required alignment positions in homology groups
    :param num_members: int - value minimal required number of members in homology groups
    :param num_uni_mem: int - value minimal required number of unique members in homology groups
    :return: list of remaining homology groups, dataframe containing alignment stats
    """
    # create dataframe of the combined align info
    df_align_info = pd.DataFrame(combined_align_info, columns=['id', 'numb_members', 'num_genomes', 'align_length',
                                                               'num_var_pos', 'num_inf_pos'])
    logging.info("Number of Homology groups before applying the config filters: {}".format(len(df_align_info)))
    # filters are done individually as to keep it readable
    # filter out short align lengths
    df_align_processed = df_align_info.query('align_length >= @align_len')
    # filter out any groups with lower number of members (number of sequences in a group)
    df_align_processed = df_align_processed.query('numb_members >= @num_members')
    # filter out groups with lower number of unique members (unique genomes)
    df_align_processed = df_align_processed.query('num_genomes >= @num_uni_mem')

    grp_ids_filtered = list(df_align_processed['id'])

    # report the number of groups back to the log file
    logging.info("Number of groups after alignment based filters from the config: {}".format(len(grp_ids_filtered)))
    if len(grp_ids_filtered) < 1:
        raise ValueError("Number of groups fell below 1 after filtering, relax the filters to resolve the issue.")

    # make list of full paths for each remaining group
    grp_path_list = []
    for grp in grp_ids_filtered:
        grp_path_list.append(os.path.join(pan_grp_path, grp))

    # make individual folders in PanVA dir for each of the remaining groups
    try:
        os.mkdir(panva_path)
    except FileExistsError:
        pass

    for grp in grp_ids_filtered:
        try:
            # join the panva_path and group
            single_panva_in_path = os.path.join(panva_path, grp)
            # Create target Directory
            os.mkdir(single_panva_in_path)
        except FileExistsError:
            pass
    # return the df (basis for homologies.json), the list of full paths to groups that passed filters
    return df_align_processed, grp_path_list
