#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:

# Imports:
import os
import pandas as pd
import numpy as np
# Functions:


def genome_annotation(single_path, panva_p):
    """
    If the reference genomes have VCF files collect this information

    :param single_path: str - path to single homology grp id in alignments/../.. .
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :return: list and Dataframe with reference genome information on CDS/exon/intron.
    """

    # get hom id and make output path
    hom_id = single_path.rsplit('/', 1)[1]
    # panva out path
    panva_out = os.path.join(panva_p, hom_id)
    # get the hom group id and add required file to path
    hom_trimm_inf_path = os.path.join(single_path, "input/trimmed.info")

    with open(hom_trimm_inf_path, 'r') as trim:
        trim_lines = trim.readlines()
    trim.close()
    # remove heading
    counter = 1
    for line in trim_lines:
        counter += 1
        # different PanTools versions pangenomes going around no consistent end nor parse line 
        if line.startswith("Total") or line.startswith("#Total"):
            break

    no_head_trim = trim_lines[counter:]

    accessions = []
    trim_start = []
    trim_stop = []
    # for i in no_head_trim:
    #     accessions.append(i.split(':')[0])
    #
    #     nrs = i.split(':')[1].strip('\n').split(',')
    #     start = nrs[0].strip(' ')
    #     end = nrs[1].strip(' ')
    #     trim_start.append(start)
    #     trim_stop.append(end)
    for i in no_head_trim:

        segments = i.strip('\n').split(',')
        acces = segments[0]
        start = segments[1]
        end = segments[2]
        accessions.append(acces)
        trim_start.append(start)
        trim_stop.append(end)


    df_trimmed = pd.DataFrame()
    df_trimmed['mRNA_id'] = accessions
    df_trimmed['start'] = trim_start
    df_trimmed['stop'] = trim_stop
    # get just the ref genomes
    lst_refs = [i for i in df_trimmed['mRNA_id'].to_list() if '|' not in i]
    # lst_refs_new = [i for i in df_trimmed['mRNA_id'].to_list() if '|' not in i]

    # subset information from df
    df_ref_trimmed = df_trimmed[df_trimmed['mRNA_id'].isin(lst_refs)]

    return lst_refs, df_ref_trimmed


def genome_annot_al_pos(single_path, panva_p, lst_refs, df_ref_trimmed):
    """
    If the reference genomes have VCF files add this level of annotation to the alignments

    :param single_path: str - path to single homology grp id in alignments/../.. .
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :param lst_refs: list - list of the reference genomes present in the homology group
    :param df_ref_trimmed: Dataframe - Dataframe with the information of the reference genome annotation
    :return:
    """

    # get hom id and make output path
    hom_id = single_path.rsplit('/', 1)[1]
    # panva out path
    panva_out = os.path.join(panva_p, hom_id)
    # get the hom group id and add required file to path
    hom_nuc_struct_path = os.path.join(single_path, "input/nuc.structure.tsv")
    if not os.path.isfile(hom_nuc_struct_path):
        hom_nuc_struct_path = os.path.join(single_path, "input/var.structure.tsv")
    # open the nuc_struct
    nuc_struct_df = pd.read_csv(hom_nuc_struct_path, delimiter='\t')

    # just open the al pos for now later fix that it takes it from memory
    all_pos = pd.read_csv(os.path.join(panva_out, 'alignments.csv'), low_memory=False)

    ref_holding_list = []
    for ref in lst_refs:
        num_o_gaps_1 = all_pos[(all_pos['mRNA_id'] == ref) & (all_pos['nucleotide'] == '-')]
        nuc_struc_1 = nuc_struct_df[nuc_struct_df['mRNA'] == ref]
        # get start
        nr_start = int(df_ref_trimmed[df_ref_trimmed['mRNA_id'] == ref]['start'])
        nr_stop = int(df_ref_trimmed[df_ref_trimmed['mRNA_id'] == ref]['stop'])

        nuc_struc_1 = nuc_struc_1.sort_values('Start').reset_index()

        nuc_struc_1['trimmed_start'] = [(nuc_struc_1['Start'][item] - nuc_struc_1['Start'][0])+1-nr_stop for item in
                                        range(len(nuc_struc_1))]
        nuc_struc_1['length'] = [(nuc_struc_1['Stop'][item] - nuc_struc_1['Start'][item]) for item in range(len(nuc_struc_1))]
        nuc_struc_1['trimmed_stop'] = [(nuc_struc_1['trimmed_start'][item] + nuc_struc_1['length'][item]) for item in range(len(nuc_struc_1))]
        # get the features

        # first version just do CDS later add others
        nuc_struc_feat_1 = nuc_struc_1[nuc_struc_1['Feature'].isin(['CDS', 'exon'])].drop_duplicates()
        # nuc_struc_exon_1 = nuc_struc_1[nuc_struc_1['Feature'].isin(['exon'])].drop_duplicates()
        lst_gaps_1 = num_o_gaps_1['position'].to_list()

        # set correct index
        nuc_struc_feat_1 = nuc_struc_feat_1.reset_index()
        nuc_struc_feat_1_sorted = nuc_struc_feat_1.sort_values('trimmed_start')

        nuc_struc_feat_1_sorted['trim_new_start'] = nuc_struc_feat_1_sorted['trimmed_start']
        nuc_struc_feat_1_sorted['trim_new_stop'] = nuc_struc_feat_1_sorted['trimmed_stop']

        df_feats = nuc_struc_feat_1_sorted
        # can be improved in speed but for now fast enough
        for entry in range(len(df_feats)):

            rnge_len = range(df_feats['trim_new_start'][entry], df_feats['trim_new_stop'][entry])
            if len(rnge_len) > 0:
                for gap_pos in lst_gaps_1:

                    if gap_pos in rnge_len or gap_pos < rnge_len[0]:

                        if gap_pos > rnge_len[0]:
                            # just end change ?
                            old_end = df_feats.at[entry, 'trim_new_stop']
                            df_feats.at[entry, 'trim_new_stop'] = old_end + 1
                        else:
                            # both end and start change
                            old_start = df_feats.at[entry, 'trim_new_start']
                            df_feats.at[entry, 'trim_new_start'] = old_start + 1

                            old_end = df_feats.at[entry, 'trim_new_stop']
                            df_feats.at[entry, 'trim_new_stop'] = old_end + 1
                    elif gap_pos > rnge_len[-1]:
                        old_end = df_feats.at[entry, 'trim_new_stop']
                        df_feats.at[entry, 'trim_new_stop'] = old_end + 1
                    else:
                        print("homid {} and gap still outside but not larger".format(hom_id, gap_pos))
                        # print(rnge_len)

        df_feats = df_feats.sort_values(by=['Feature'])

        ref_pos_cds_list = []
        ref_pos_exon_list = []
        ref_feat_exon_list = []
        ref_feat_cds_list = []

        for entry in range(len(df_feats)):
            if df_feats['Feature'][entry] == 'CDS':
                cds_range_pos = [*range(df_feats['trim_new_start'][entry], df_feats['trim_new_stop'][entry]+1)]
                ref_pos_cds_list.extend(cds_range_pos)
                ref_feat_cds_list.extend([True]*len(cds_range_pos))

            elif df_feats['Feature'][entry] == 'exon':
                exon_range_pos = [*range(df_feats['trim_new_start'][entry], df_feats['trim_new_stop'][entry] + 1)]
                ref_pos_exon_list.extend(exon_range_pos)
                ref_feat_exon_list.extend([True]*len(exon_range_pos))

        cds_data = {'position': ref_pos_cds_list, 'cds': ref_feat_cds_list}
        new_idea_cds_df = pd.DataFrame(cds_data)
        new_idea_cds_df = new_idea_cds_df.drop_duplicates()

        exon_data = {'position': ref_pos_exon_list, 'exon': ref_feat_exon_list}
        new_idea_exon_df = pd.DataFrame(exon_data)
        new_idea_exon_df = new_idea_exon_df.drop_duplicates()

        new_idea_df = pd.merge(new_idea_cds_df, new_idea_exon_df, on='position', how='outer')
        # new_idea_df = pd.merge(ref_pos_cds_df, ref_pos_exon_df, on='position', how='outer')
        new_idea_df['position'] = pd.to_numeric(new_idea_df['position'])
        new_idea_df = new_idea_df.drop_duplicates()

        # FILL intron positions
        end_of_pos = max(new_idea_df['position'])
        full_range = [*range(1, end_of_pos+1)]

        new_idea_df = new_idea_df.set_index('position').reindex(full_range).fillna(False).reset_index()

        new_idea_df['mRNA_id'] = ref

        new_idea_df = new_idea_df.sort_values(by=['position'])

        max_nr = max(list(set(nuc_struc_1['length'].to_list())))
        length_aligned = max_nr - nr_start - nr_stop + len(lst_gaps_1) + 1
        # print('length_aligned', length_aligned)
        df_coding_1_trimmed = new_idea_df[new_idea_df['position'].isin(range(1, length_aligned + 1))]

        ref_holding_list.append(df_coding_1_trimmed)

    df_ref_feats = pd.concat(ref_holding_list)

    df_all_ref_pos = all_pos[all_pos['mRNA_id'].isin(lst_refs)]

    df_nuc_ref_struct = df_ref_feats.merge(df_all_ref_pos, on=['mRNA_id', 'position'])
    df_nuc_ref_struct = df_nuc_ref_struct.drop_duplicates()
    df_nuc_ref_struct.to_csv(os.path.join(panva_out, "annotations.csv"), index=False)

    return


def pool3wrap(single_id_path, panva_path):
    """
    Wrapper function to make use of multiprocessing. executes genome_annotation.

    :param single_id_path:
    :param panva_path:
    :return:
    """

    ref_list, df_refs = genome_annotation(single_id_path, panva_path)
    genome_annot_al_pos(single_id_path, panva_path, ref_list, df_refs)

    return
