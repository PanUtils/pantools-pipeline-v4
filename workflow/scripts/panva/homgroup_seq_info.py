#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:
# For each homology group keep track of sequence_id, genome_nr, functioning as a link.
# keeping all essential information together.

# Imports:
import pandas as pd
import os
from functools import reduce


# Functions:
def create_seq_info(single_id_path, panva_p):
    """
    Collects information on the individual sequences in a homology group.

    :param single_id_path: str - path to single homology grp id in alignments/../.. .
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :return: df_seq_info - dataframe - contains information on sequences in the homology group
    """
    # get single id code
    single_id = single_id_path.rsplit('/', 1)[1]
    # add extension to path to required file
    single_id_seq_inf = os.path.join(single_id_path, "input/sequences.info")

    with open(single_id_seq_inf, 'r') as rd:
        comment_lines = rd.readlines()
    rd.close()

    # finding the headers
    index_log = 0
    for line in comment_lines:
        index_log = index_log + 1
        if line.startswith('#genome'):
            break

    line_info = comment_lines[index_log:]
    lst_sequences = []
    # sequence information section
    for string in line_info:
        string = string.strip('\n').split(', ')
        # print('line: ', string)
        lst_sequences.append(string)

    df_seq_info = pd.DataFrame(lst_sequences, columns=['genome_nr', 'gene_name', 'mRNA_id', 'node_id', 'address',
                                                       'strand'])

    df_seq_info.drop(0)

    # Check for accessions
    single_id_acc_inf = os.path.join(single_id_path, "input/genome_order.info")
    # naming column to mRNA_id for merging purposes in the even there is resequenced accessions in the data
    seq_acc_info = pd.read_csv(single_id_acc_inf, header=0, names=['genome_nr', 'mRNA_id'])

    df_seq_info = pd.merge(df_seq_info, seq_acc_info, on=['genome_nr', 'mRNA_id'], how='outer')
    df_seq_info.to_csv(os.path.join(os.path.join(panva_p, single_id), "sequence.info"), index=False)

    return df_seq_info


def seq_df_maker(fasta_path):
    """
    Collects the fastas from the homology groups msa.

    :param fasta_path: str - path to the fasta file
    :return: data entries of the fasta file
    """

    with open(fasta_path, 'r') as f:
        nuc_lines = f.readlines()
    f.close()

    nuc_lines = [i.strip('\n') for i in nuc_lines]
    new_lines = list(enumerate(nuc_lines))

    # get header lines
    headers = []
    [headers.append(new_lines[i]) for i in range(len(new_lines)) if new_lines[i][1].startswith('>')]

    seqs = []
    # calculate nr of lines between headers fasta
    nr_lines = headers[1][0] - headers[0][0]

    for i in range(0, len(headers)):
        idx = headers[i][0]
        # print(idx, idx+nr_lines)

        line_seq = [lin[1] for lin in new_lines[idx + 1:idx + nr_lines]]
        seqs.append(''.join(line_seq))
    # align_length = len(seqs)
    headers_clean = []

    # omitting the address because already in sequences json
    for header in headers:
        headers_clean.append(header[1].split(' ')[0].split('>')[1])

    # create dataframe with header and sequences (nuc trimmed)
    df_sequences = pd.DataFrame()
    df_sequences['mRNA_id'] = headers_clean
    sequence_type = fasta_path.rsplit('.', 1)[0]
    sequence_type = sequence_type.rsplit('/', 1)[1]

    df_sequences['{}_seq'.format(sequence_type)] = seqs

    return df_sequences


def merge_sequences(single_id_path, panva_path, method='nuc_trimmed'):
    """
    Combines all sequence information in a homology group into a single dataframe.

    :param single_id_path: str - path to single homology grp id in alignments/../.. .
    :param panva_path: str - full path to output dir (serves as input for PanVA).
    :param method: str - define type of sequences (default='nuc_trimmed')
    :return: dataframe of all sequences in the homology group
    """

    # select individual id
    single_id = single_id_path.rsplit('/', 1)[1]
    # direct output to panva
    panva_out = os.path.join(panva_path, single_id)

    allowed_fasta = ["nuc.fasta", "nuc_trimmed.fasta", "prot.fasta", "prot_trimmed.fasta", "var.fasta", "var_trimmed.fasta"]
    single_id_output_path = os.path.join(single_id_path, "output/")

    hold_fasta_path = []
    if method == 'all':
        for file in allowed_fasta:
            if os.path.isfile(os.path.join(single_id_output_path, file)):
                hold_fasta_path.append(os.path.join(single_id_output_path, file))
            else:
                pass
    elif method == 'nuc_trimmed':
        hold_fasta_path.append(os.path.join(single_id_path, "output/nuc_trimmed.fasta"))
    elif method == 'var_trimmed':
        hold_fasta_path.append(os.path.join(single_id_path, "output/var_trimmed.fasta"))
    else:
        raise ValueError("The method for fasta sequence merging should be 'all' or 'nuc_trimmed'.")

    all_fasta_list = []
    for file in hold_fasta_path:
        tmp_var = seq_df_maker(file)
        all_fasta_list.append(tmp_var)

    df_all_sequences = reduce(lambda left, right: pd.merge(left, right, on='mRNA_id', how='outer'), all_fasta_list)

    # with the newest version of PanTools new prefixes are introduced however PanVA does not have these
    # Therefore here we have to overwrite the sequencetype (for now)
    if method == 'var_trimmed':
        df_all_sequences.rename(columns={'var_trimmed_seq': 'nuc_trimmed_seq'}, inplace=True)

    df_all_sequences.to_csv(os.path.join(panva_out, 'sequences.csv'), index=False)

    return df_all_sequences
