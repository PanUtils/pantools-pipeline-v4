"""
Run PanTools subcommands for phased pangenomics.
Contains the following subcommands:

add_phasing
add_repeats
"""

rule all_phasing:
    input:
        f"{config['results']}/done/pangenome.add_phasing.done",
        f"{config['results']}/done/pangenome.add_repeats.done"

rule add_phasing:
    """Add phasing information to the pangenome."""
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
        workflow.cores * 0.6
    shell:
        "{pantools} add_phasing {params.opts} {params.database} {input.phasing}"

rule add_repeats:
    """Add repeats to the pangenome."""
    input:
        "{results}/done/validate.annotations.done",
        "{results}/done/pangenome.build_pangenome.done",
        annotations = config['annotations'],
    output:
        touch("{results}/done/pangenome.add_repeats.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['add_repeats.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.add_repeats.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} add_repeats -f {params.opts} {params.database} {input.annotations}"

