#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:
# Collect the data from the PanTools pangenome based on the config file requirements.

# Imports:
import logging
import os

# Function(s):


def collect_data(pan_grp_path, selection):
    """
    List's all homologies based on selection checking if the group possess all required files.

    :param pan_grp_path: str - path to the PanTools database alignment/<msa_type>/grouping_v*/
    :param selection: either path to file.txt with homology group id's (sep=,) or 'all'
    :return: list of the full paths of homology groups that made it through initial filters
    """

    if isinstance(selection, str):
        if selection == 'all':
            # keep information for the log file
            logging.info("Selection method for groups is 'all' from: {}".format(pan_grp_path))
            # create list of all selected groups
            grp_ids = [folder for folder in os.listdir(pan_grp_path) if os.path.isdir(os.path.join(pan_grp_path, folder))]
        elif os.path.isfile(selection):
            logging.info("Selection method for groups is 'file': {}".format(pan_grp_path))
            grp_ids = []
            with open(selection, 'r') as grp_file:
                for line in grp_file:
                    if ',' in line:
                        hold_ids = [elem.strip() for elem in line.split(',')]
                        grp_ids.extend(hold_ids)
                    else:
                        grp_ids = [line.strip()for line in grp_file]
            grp_file.close()
        else:
            logging.error("The selection as given in the config file should be 'all' or an 'full/path/'.")
            raise ValueError("The selection in the config should be 'all' or 'full/path'.")

    # cross-check the selection for variant positions.
    hom_w_var_pos = os.path.join(pan_grp_path, 'groups_with_var_inf_positions_nuc.txt')
    hom_w_exc_trim = os.path.join(pan_grp_path, 'groups_excluded_based_on_trimming.txt')

    # number of groups before exclusion based on trimming.
    pre_excl_num = len(grp_ids)
    logging.info("Number of groups selected before check of exclusion based on trimming is: {}".format(pre_excl_num))
    # check if any homology id groups are excluded based on trimming.
    if os.path.exists(hom_w_exc_trim):
        with open(hom_w_exc_trim, 'r') as ex_t:
            for line in ex_t:
                remove_id = [entry.strip() for entry in line.split(',')]
        ex_t.close()
        grp_ids_clean = [x for x in grp_ids if x not in remove_id]
        # rename to original list
        grp_ids = grp_ids_clean

    # Number of groups after exclusion based on trimming
    pre_filt_num = len(grp_ids)
    logging.info("Number of groups after check of exclusion based on trimming is: {}".format(pre_filt_num))

    # Check which homology group id don't have var inf positions
    if os.path.exists(hom_w_var_pos):
        with open(hom_w_var_pos, 'r') as var_p:
            all_lines = var_p.readlines()
        var_p.close()
        # select second line as this contains the homology ids with var positions
        second_line = all_lines[1]
        keep_ids = [entry.strip() for entry in second_line.split(',')]
        grp_ids_var_clean = [x for x in grp_ids if x in keep_ids]
        # back to regular list naming convention
        grp_ids = grp_ids_var_clean

    # Number of groups after variant aligned positions check
    post_filt_num = len(grp_ids)
    logging.info("Number of groups after check of presence of variant positions: {}".format(post_filt_num))
    if post_filt_num < 1:
        logging.error("The remaining number of groups fell below 0.")
        raise ValueError("Can't pre-process less than 1 group, please increase the selection or alter filtering")
    # Combine the remaining group ids with the full path of the pangenome database. stored in list
    grp_path_list = []
    for group in grp_ids:
        grp_path_list.append(os.path.join(pan_grp_path, group))
    # list containing the remaining groups after initial filters.
    return grp_path_list
