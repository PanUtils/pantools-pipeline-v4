"""
Run PanTools subcommands for classification.
Contains the followig subcommands:

functional_classification
function_overview
gene_classification
go_enrichment
grouping_overview
kmer_classification
metrics
pangenome_structure  
"""

rule functional_classification:
    """Classify functional annotations as core, accessory or unique."""
    input:
        "{results}/done/{type}.add_functions.done",
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
    output:
        done = touch("{results}/done/{type}.functional_classification.done"),
        file1 = "{results}/{type}_db/function/functional_classification/functional_classification_overview.txt",
        file2 = "{results}/{type}_db/function/functional_classification/core_functions.txt",
        file3 = "{results}/{type}_db/function/functional_classification/accessory_functions.txt",
        file4 = "{results}/{type}_db/function/functional_classification/unique_functions.txt",
    params:
        database = "{results}/{type}_db",
        opts = config['functional_classification.opts'],
    benchmark:
        "{results}/benchmarks/{type}.functional_classification.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} functional_classification {params.opts} {params.database}"

rule function_overview:
    """Create an overview table for each functional annotation type in the pangenome."""
    input:
        "{results}/done/{type}.add_functions.done"
    output:
        done = touch("{results}/done/{type}.function_overview.done"),
        file1 = "{results}/{type}_db/function/functions_per_group_and_mrna.csv",
        file2 = "{results}/{type}_db/function/function_counts_per_group.csv",
    params:
        database = "{results}/{type}_db",
        opts = config['function_overview.opts'],
    benchmark:
        "{results}/benchmarks/{type}.function_overview.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} function_overview -f {params.opts} {params.database}"

rule gene_classification:
    """Classify the gene repertoire as core, accessory or unique."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/done/{type}.add_pavs.done" if config['pav'] else [],
    output:
        touch("{results}/done/{type}.gene_classification.done"),
        "{results}/{type}_db/gene_classification/accessory_combinations.csv",
        "{results}/{type}_db/gene_classification/accessory_groups.csv",
        "{results}/{type}_db/gene_classification/all_homology_groups.csv",
        "{results}/{type}_db/gene_classification/classified_groups.csv",
        "{results}/{type}_db/gene_classification/cnv_core_accessory.txt",
        "{results}/{type}_db/gene_classification/core_groups.csv",
        "{results}/{type}_db/gene_classification/gene_classification_overview.txt",
        "{results}/{type}_db/gene_classification/gene_distance.tree",
        "{results}/{type}_db/gene_classification/group_size_occurrences.txt",
        "{results}/{type}_db/gene_classification/shared_unshared_gene_count.csv",
        "{results}/{type}_db/gene_classification/single_copy_orthologs.csv",
        "{results}/{type}_db/gene_classification/unique_groups.csv",
    params:
        database = "{results}/{type}_db",
        opts = config['gene_classification.opts'],
    benchmark:
        "{results}/benchmarks/{type}.gene_classification.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        """
        {pantools} gene_classification -f {params.opts} {params.database}
        Rscript {params.database}/gene_classification/gene_distance_tree.R
        """

rule find_dispensable_homology_groups:
    """Find dispensible homology groups to use for go_enrichment."""
    input:
        accessory = "{results}/{type}_db/gene_classification/accessory_groups.csv",
        unique = "{results}/{type}_db/gene_classification/unique_groups.csv",
    output:
        "{results}/{type}_db/gene_classification/dispensable_groups.csv",
    shell:
        "awk 'BEGIN{{ORS = \",\";}} $1 !~ /^[G#]/ && $1 !~ /^$/' {input} | sed 's/,$/\\n/g' > {output}"

rule go_enrichment:
    """Identify over- or underrepresented GO terms in a set of genes."""
    input:
        "{results}/done/{type}.add_functions.done",
        groups = "{results}/{type}_db/gene_classification/dispensable_groups.csv",
    output:
        touch("{results}/done/{type}.go_enrichment.done")
    params:
        database = "{results}/{type}_db",
        opts = config['go_enrichment.opts'],
    benchmark:
        "{results}/{type}_db/benchmarks/{type}.go_enrichment.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} go_enrichment -f -H={input.groups} {params.opts} {params.database}"

rule grouping_overview:
    """Create an overview table for every homology grouping in the pangenome."""
    input:
        "{results}/done/{type}.grouping.done"
    output:
        touch("{results}/done/{type}.grouping_overview.done"),
        "{results}/{type}_db/group/grouping_overview.txt",
    params:
        database = "{results}/{type}_db",
        opts = config['grouping_overview.opts'],
    benchmark:
        "{results}/benchmarks/{type}.grouping_overview.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} grouping_overview {params.opts} {params.database}"

rule kmer_classification:
    """Calculate the number of core, accessory, unique k-mer sequences."""
    input:
        "{results}/done/pangenome.build_pangenome.done"
    output:
        touch("{results}/done/pangenome.kmer_classification.done"),
        "{results}/pangenome_db/kmer_classification/genome_kmer_distance.tree",
        "{results}/pangenome_db/kmer_classification/kmer_classification_overview.txt",
    params:
        database = "{results}/pangenome_db",
        opts = config['kmer_classification.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.kmer_classification.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        """
        {pantools} kmer_classification -f {params.opts} {params.database}
        Rscript {params.database}/kmer_classification/genome_kmer_distance_tree.R
        """

rule metrics:
    """Generates relevant metrics of the pangenome."""
    input:
        "{results}/done/{type}.build_{type}.done"
    output:
        touch("{results}/done/{type}.metrics.done"),
        "{results}/{type}_db/metrics/metrics.txt",
        "{results}/{type}_db/metrics/metrics_per_genome.csv",
    params:
        database = "{results}/{type}_db",
        opts = config['metrics.opts'],
    benchmark:
        "{results}/benchmarks/{type}.metrics.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        "{pantools} metrics -f {params.opts} {params.database}"
    
rule pangenome_structure:
    """Determine the openness of the pangenome based on homology groups."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_pavs.done" if config['pav'] else [],
    output:
        touch("{results}/done/{type}.pangenome_structure.done"),
        "{results}/{type}_db/pangenome_size/gene/pangenome_size.txt",
    params:
        database = "{results}/{type}_db",
        opts = config['pangenome_structure.opts'],
    benchmark:
        "{results}/benchmarks/{type}.pangenome_structure.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        {pantools} pangenome_structure -f --threads={threads} {params.opts} {params.database}
        for file in {params.database}/pangenome_size/gene/*.R; do Rscript $file; done
        """

rule kmer_structure:
    """Determine the openness of the pangenome based on k-mer sequences."""
    input:
        "{results}/done/pangenome.grouping.done",
    output:
        touch("{results}/done/pangenome.pangenome_structure_kmer.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['kmer_structure.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.pangenome_structure_kmer.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        {pantools} pangenome_structure --kmer --threads={threads} {params.opts} {params.database}
        for file in {params.database}/pangenome_size/kmer/*.R; do Rscript $file; done
        """
