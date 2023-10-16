"""
Prepare PanTools output for PanVa.
"""

rule to_panva:
    input:
        lambda wildcards: msa_done(wildcards.type),
        "{results}/done/{type}.ani.done",
        "{results}/done/{type}.gene_classification.done",
        "{results}/done/{type}.group_info.done",
        "{results}/done/{type}.kmer_classification.done"
    output:
        touch("{results}/done/{type}.panva.done")
    params:
        panva_config = config['panva_config']
    conda:
        "../envs/pantova.yaml"
    shell:
        "python3 workflow/scripts/panva/pan_to_va.py {params.panva_config}"

rule panva:
    """""Run all PanTools functions necessary to create a PanVa instance."""
    input:
        "{}/done/pangenome.panva.done".format(config['results'])