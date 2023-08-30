"""
Prepare PanTools output for PanVa.
"""

rule panva:
    """""Run all PanTools functions necessary to create a PanVa instance."""
    input:
        "{}/done/pangenome.ani.done".format(config['results']),
        "{}/done/pangenome.gene_classification.done".format(config['results']),
        "{}/done/pangenome.group_info.done".format(config['results']),
        "{}/done/pangenome.kmer_classification.done".format(config['results']),
        "{}/done/pangenome.msa.done".format(config['results'])