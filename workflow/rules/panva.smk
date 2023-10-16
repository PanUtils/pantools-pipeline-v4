"""
Prepare PanTools output for PanVa.
"""

rule to_panva:
    input:
        lambda wildcards: msa_done("pangenome"),
        "{results}/done/pangenome.ani.done",
        "{results}/done/pangenome.gene_classification.done",
        "{results}/done/pangenome.group_info.done",
        "{results}/done/pangenome.kmer_classification.done"
    output:
        touch("{results}/done/panva.done")
    params:
        panva_config = config['panva_config']
    conda:
        "../envs/pantova.yaml"
    shell:
        "python3 ../scripts/pant_to_va.py {params.panva_config}"

rule panva:
    """""Run all PanTools functions necessary to create a PanVa instance."""
    input:
        "{}/done/panva.done".format(config['results'])