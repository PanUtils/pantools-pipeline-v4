"""
Run PanTools subcommands for construction.
Contains the following subcommands:

add_annotations
add_functions
add_phenotypes
build_pangenome
build_panproteome
busco_protein
change_grouping
group
optimal_grouping 
"""

rule add_annotations:
    """Add annotations to the pangenome."""
    input:
        "{results}/done/validate.annotations.done",
        "{results}/done/pangenome.build_pangenome.done",
        annotations = config['annotations'],
    output:
        touch("{results}/done/pangenome.add_annotations.done")
    params:
        database = f"{config['construction']}/pangenome_db",
        opts = config['add_annotations.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.add_annotations.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} add_annotations -f {params.opts} {params.database} {input.annotations}"
                
rule add_functions:
    """"Add functional annotations to the pangenome."""
    input:
        "{results}/done/validate.functions.done",
        lambda wildcards: proteins_done(wildcards.type),
        functions = config['functions']
    output:
        touch("{results}/done/{type}.add_functions.done")
    params:
        database = f"{config['construction']}/{{type}}_db",
        opts = config['add_functions.opts'],
    benchmark:
        "{results}/benchmarks/{type}.add_functions.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} add_functions -f -F=resources/functional_databases {params.opts} {params.database} {input.functions}"
        
rule add_phenotypes:
    """Add phenotype data to the pangenome."""
    input:
        "{results}/done/{type}.build_{type}.done",
        phenotypes = config['phenotypes']
    output:
        touch("{results}/done/{type}.add_phenotypes.done"),
    params:
        database = f"{config['construction']}/{{type}}_db",
        opts = config['add_phenotypes.opts'],
    benchmark:
        "{results}/benchmarks/{type}.add_phenotypes.txt"
    conda:
        "../envs/pantools.yaml",
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} add_phenotypes {params.opts} {params.database} {input.phenotypes}"

rule build_pangenome:
    """Build a pangenome from a set of genomes."""
    input:
        "{results}/done/validate.genomes.done",
        genomes = config['genomes']
    output:
        touch("{results}/done/pangenome.build_pangenome.done"),
    params:
        database = f"{config['construction']}/pangenome_db",
        opts = config['build_pangenome.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.build_pangenome.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} build_pangenome -f -t={threads} --scratch-directory={temp_dir} {params.opts} {params.database} {input.genomes}"

rule build_panproteome:
    """Build a panproteome from a set of proteins."""
    input:
        "{results}/done/validate.proteomes.done",
        proteomes = config['proteomes']
    output:
        touch("{results}/done/panproteome.build_panproteome.done"),
    params:
        database = f"{config['construction']}/panproteome_db",
        opts = config['build_panproteome.opts'],
    benchmark:
        "{results}/benchmarks/panproteome.build_panproteome.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} build_panproteome -f {params.opts} {params.database} {input.proteomes}"

rule busco_protein:
    """Identify BUSCO genes in the pangenome."""
    input:
        lambda wildcards: proteins_done(wildcards.type)
    output:
        touch("{results}/done/{type}.busco_protein.{busco}.done")
    params:
        database = f"{config['construction']}/{{type}}_db",
        opts = config['busco_protein.opts'],
    benchmark:
        "{results}/benchmarks/{type}.busco_protein.{busco}.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} busco_protein -t={threads} {params.opts} --busco10={wildcards.busco} {params.database}"

rule change_grouping:
    """"Change the active grouping based on the result of optimal_grouping."""
    input:
        "{results}/done/{type}.optimal_grouping.done"
    output:
        touch("{results}/done/{type}.change_grouping.done")
    params:
        database = f"{config['construction']}/{{type}}_db",
        opts = config['add_functions.opts'],
    benchmark:
        "{results}/benchmarks/{type}.change_grouping.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        """
        relaxation=$(sort -k14 -n {params.database}/optimal_grouping/grouping_overview.csv | head -1 -c 2 | tail -c 1)
        {pantools} change_grouping -v=$relaxation {params.opts} {params.database}
        """
        
rule group:
    """Generate homology groups based on similarity of protein sequences."""
    input:
        lambda wildcards: proteins_done(wildcards.type)
    output:
        touch("{results}/done/{type}.group.done")
    params:
        database = f"{config['construction']}/{{type}}_db",
        relaxation = config["group.relaxation"],
        opts = config['group.opts'],
    benchmark:
        "{results}/benchmarks/{type}.group.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} group -f -t={threads} {params.opts} --relaxation={params.relaxation} {params.database}"

rule optimal_grouping:
    """Find the most suitable settings for group using BUSCO output."""
    input:
        "{{results}}/done/{{type}}.busco_protein.{busco}.done".format(busco=config['busco_protein.odb10'])
    output:
        touch("{results}/done/{type}.optimal_grouping.done")
    params:
        database = f"{config['construction']}/{{type}}_db",
        busco = config["busco_protein.odb10"],
        opts = config['optimal_grouping.opts'],
    benchmark:
        "{results}/benchmarks/{type}.optimal_grouping.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        {pantools} optimal_grouping -t={threads} {params.opts} {params.database} {params.database}/busco/{params.busco}
        Rscript {params.database}/optimal_grouping/optimal_grouping.R
        """
        
checkpoint grouping:
    """Checkpoint for an active grouping either using a given relaxation or BUSCO set"""
    input:
        "{results}/done/{type}.group.done" if config['group.relaxation'] else \
            "{results}/done/{type}.change_grouping.done"
    output:
        touch("{results}/done/{type}.grouping.done")

rule get_group_ids:
    input:
        "{results}/done/{type}.grouping.done",
        gene_selection = config["gene_selection"]
    output:
        "{results}/{type}_db/homology_selection.txt"
    params:
        all_groups = f"{config['construction']}/{{type}}_db/pantools_homology_groups.txt",
    script:
        "../scripts/homology_selection.sh"

checkpoint construction:
    """Checkpoint for all pangenome construction."""
    input:
        lambda wildcards: construction_done(wildcards.type)
    output:
        touch("{results}/done/{type}.construction.done")
    params:
        database = f"{config['construction']}/{{type}}_db",
        results = config['results']
    shell:
        "mv {params.database} {params.results}"