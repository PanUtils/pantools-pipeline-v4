"""
Run PanTools subcommands for phased pangenomics.
Contains the following subcommands:

add_phasing
add_repeats
repeat_overview
calculate_synteny
add_synteny
synteny_statistics
"""

rule all_phasing:
    input:
        f"{config['results']}/done/pangenome.add_phasing.done",
        f"{config['results']}/done/pangenome.repeat_overview.done",
        f"{config['results']}/done/pangenome.synteny_overview.done",
        f"{config['results']}/done/pangenome.blast.done"

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

rule repeat_overview:
    """Provide repeat statistics."""
    input:
        "{results}/done/pangenome.add_repeats.done"
    output:
        touch("{results}/done/pangenome.repeat_overview.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['repeat_overview.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.repeat_overview.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} repeat_overview -f {params.opts} {params.database}"

rule calculate_synteny:
    """Calculate synteny information using MCSCanX."""
    input:
        "{results}/done/pangenome.group.done",
    output:
        touch("{results}/done/pangenome.calculate_synteny.done"),
        synteny = "{results}/pangenome_db/synteny/mcscanx.collinearity"
    params:
        database = "{results}/pangenome_db",
        opts = config['calculate_synteny.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.calculate_synteny.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} calculate_synteny -f {params.opts} -t {threads} {params.database}"

rule add_synteny:
    """Add synteny information to the pangenome."""
    input:
        "{results}/done/pangenome.group.done",
        synteny = config['syneny'] if config['synteny'] else "{results}/pangenome_db/synteny/mcscanx.collinearity"
    output:
        touch("{results}/done/pangenome.add_synteny.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['add_synteny.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.add_synteny.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} add_synteny -f {params.opts} {params.database} {input.synteny}"

rule synteny_overview:
    """Provide synteny statistics."""
    input:
        "{results}/done/pangenome.add_synteny.done",
        synteny = config['syneny'] if config['synteny'] else "{results}/pangenome_db/synteny/mcscanx.collinearity"
    output:
        touch("{results}/done/pangenome.synteny_overview.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['synteny_overview.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.synteny_overview.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} synteny_overview -f {params.opts} {params.database} {input.synteny}"

rule blast:
    """Run BLAST."""
    input:
        "{results}/done/pangenome.build_pangenome.done",
        blast = config['blast']
    output:
        touch("{results}/done/pangenome.blast.done")
    params:
        database = "{results}/pangenome_db",
        opts = config['blast.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.blast.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} blast -f {params.opts} {params.database} {input.blast}"

