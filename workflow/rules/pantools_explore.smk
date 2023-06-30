"""
Run PanTools subcommands for pangenome exploration.
Contains the followig subcommands:

compare_go
group_info
show_go
"""

rule show_go:
    """
    For a given GO term, show the child terms, all parent terms 
    higher in the hierarchy, and connected mRNA nodes.
    """
    input:
        "{results}/done/{type}.add_functions.done"
    output:
        done = touch("{results}/done/{type}.show_go.done")
    params:
        database = "{results}/{type}_db",
        opts = config['show_go.opts'],
    benchmark:
        "{results}/benchmarks/{type}.show_go.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} show_go {params.opts}{params.database}"

rule compare_go:
    """For two given GO terms, move up in the GO hierarchy to see if they are related."""
    input:
        "{results}/done/{type}.add_functions.done"
    output:
        done = touch("{results}/done/{type}.compare_go.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['compare_go.opts'],
    benchmark:
        "{results}/benchmarks/{type}.compare_go.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} compare_go {params.opts} {params.database}"

rule group_info:
    """Report all available information of one or multiple homology groups."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/done/{type}.add_functions.done" if config['functions'] else [],
    output:
        done = touch("{results}/done/{type}.group_info.done"),
    params:
        database = "{results}/{type}_db",
        opts = config['group_info.opts'],
    benchmark:
        "{results}/benchmarks/{type}.compare_go.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} group_info {params.opts} {params.database}"