"""
Run PanTools subcommands for pangenome multiple sequence alignment.
Contains the followig subcommands:

msa
"""

rule msa:
    """Create multiple sequence alignments."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_functions.done" if config['functions'] else [],
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/done/{type}.add_variants.done" if config['vcf'] else [],
        "{results}/done/{type}.add_pavs.done" if config['pav'] else [],
        "{results}/{type}_db/homology_selection.txt" if config['gene_selection'] else []
    output:
        done = touch("{results}/done/{type}.msa.done"),
    params:
        database = "{results}/{type}_db",
        opts = config['msa.opts'],
        homology = "-H={results}/{type}_db/homology_selection.txt" if config['gene_selection'] else ""
    benchmark:
        "{results}/benchmarks/{type}.msa.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} msa -t={threads} {params.opts} {params.homology} {params.database}"
