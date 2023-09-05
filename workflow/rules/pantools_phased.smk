"""
Run PanTools subcommands for phased pangenomics.
Contains the following subcommands:

add_phasing
"""

rule all_phasing:
    input:
        f"{config['results']}/done/pangenome.add_phasing.done"

rule add_phasing:
    """Create multiple sequence alignments."""
    input:
        "{results}/done/pangenome.build_pangenome.done",
        phasing = config['phasing']
    output:
        done = touch("{results}/done/pangenome.add_phasing.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['add_phasing.opts'],
    benchmark:
        "{results}/benchmarks/pantools.add_phasing.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} add_phasing {params.opts} {params.database} {input.phasing}"

