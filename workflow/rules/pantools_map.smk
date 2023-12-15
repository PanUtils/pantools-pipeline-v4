"""
Run PanTools subcommands for pangenome read mapping.
Contains the following subcommands:

map
"""

rule map:
    """Map single or paired-end short reads to one or multiple genomes in the pangenome."""
    input:
        "{results}/done/pangenome.construction.done",
        sr1 = config['short_read_1'],
        sr2 = config['short_read_2']
    output:
        touch("{results}/done/pangenome.map.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['map.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.map.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} map -f -t={threads} {params.opts} {params.database} {input.sr1} {input.sr2}"