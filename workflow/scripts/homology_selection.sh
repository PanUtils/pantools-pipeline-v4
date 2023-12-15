#!/bin/bash

touch temp_homology_groups.txt
while IFS= read -r line
do
  grep "$line" ${snakemake_params[all_groups]} | sed 's/:.*$//g' >> temp_homology_groups.txt
done <  ${snakemake_input[gene_selection]}
tr '\n' ',' < temp_homology_groups.txt > ${snakemake_output[0]}
rm temp_homology_groups.txt