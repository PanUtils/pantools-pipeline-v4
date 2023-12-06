"""
Run PanTools subcommands for adding variation to the pangenome.
Contains the following subcommands:

add_pavs
add_variants
variation_overview
"""

rule add_variants:
    """Add variation from vcf files to the pangenome."""
    input:
        "{results}/done/validate.variants.done",
        "{results}/done/pangenome.add_annotations.done",
        variants = config['vcf'],
    output:
        done = touch("{results}/done/pangenome.add_variants.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['add_variants.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.add_variants.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} add_variants -t={threads} {params.opts} {params.database} {input.variants}"

rule add_pavs:
    """Add presence-absence variation (PAV) to the pangenome."""
    input:
        "{results}/done/validate.pavs.done",
        "{results}/done/{type}.add_annotations.done",
        pavs = config['pav'],
    output:
        done = touch("{results}/done/{type}.add_pavs.done")
    params:
        database = "{results}/{type}_db",
        opts = config['add_pavs.opts'],
    benchmark:
        "{results}/benchmarks/{type}.add_pavs.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} add_pavs {params.opts} {params.database} {input.pavs}"

rule variation_overview:
    """"Write an overview of all accessions (PAV and VCF) added to the pangenome."""
    input:
        "{results}/done/{type}.construction.done"
    output:
        done = touch("{results}/done/{type}.variation_overview.done")
    params:
        database = "{results}/{type}_db",
        opts = config['variation_overview.opts'],
    benchmark:
        "{results}/benchmarks/{type}.variation_overview.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} variation_overview {params.opts} {params.database}"