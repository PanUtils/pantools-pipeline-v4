"""
Run PanTools subcommands for phylogeny in the pangenome.
Contains the followig subcommands:

ani
consensus_tree
core_phylogeny
mlsa_find_genes
mlsa_concatenate
mlsa
"""

rule core_phylogeny:
    """Create a SNP tree from single-copy genes."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/done/{type}.add_variants.done" if config['vcf'] else [],
        "{results}/done/{type}.add_pavs.done" if config['pav'] else [],
        homology = "{results}/{type}_db/homology_selection.txt" if config['gene_selection'] \
            else "{results}/{type}_db/gene_classification/single_copy_orthologs.csv"
    output:
        "{results}/{type}_db/core_snp_tree/sites_per_group.csv",
        done = touch("{results}/done/{type}.core_phylogeny.done"),
    params:
        database = "{results}/{type}_db",
        mode = config['core_phylogeny.mode'],
        opts = config['core_phylogeny.opts'],
    benchmark:
        "{results}/benchmarks/{type}.core_phylogeny.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        {pantools} core_phylogeny --threads={threads} -m={params.mode} -H={input.homology} {params.opts} {params.database}
        if [ {params.mode} = NJ ]; then 
            Rscript {params.database}/core_snp_tree/core_snp_NJ_tree.R
        else
            iqtree -T AUTO -ntmax	{threads} -s {params.database}/core_snp_tree/informative.fasta -redo -bb 1000
        fi
        """

rule consensus_tree:
    """Create a consensus tree by combining gene trees from homology groups using ASTRAL-Pro."""
    input:
        "{results}/done/pangenome.gene_classification.done",
        "{results}/done/{type}.add_variants.done" if config['vcf'] else [],
        "{results}/done/{type}.add_pavs.done" if config['pav'] else [],
        "{results}/{type}_db/homology_selection.txt" if config['gene_selection'] else []
    output:
        done = touch("{results}/done/{type}.consensus_tree.done"),
    params:
        database = "{results}/{type}_db",
        opts = config['consensus_tree.opts'],
        homology = "-H={results}/{type}_db/homology_selection.txt" if config['gene_selection'] else ""
    benchmark:
        "{results}/benchmarks/{type}.consensus_tree.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} consensus_tree --threads={threads} {params.opts} {params.homology} {params.database}"

rule ani:
    """Calculate Average Nucleotide Identity (ANI) scores between genomes."""
    input:
        "{results}/done/pangenome.build_pangenome.done",
        "{results}/done/pangenome.add_phenotypes.done" if config['phenotypes'] else [],
    output:
        done = touch("{results}/done/pangenome.ani.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['ani.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.ani.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        {pantools} ani --threads={threads} {params.opts} {params.database}
        Rscript {params.database}/ANI/*/ANI_tree.R
        """
        
rule mlsa_find_genes:
    """Step 1/3 of mlsa. Search and filter suitable genes for the mlsa."""
    input:
        "{results}/done/pangenome.gene_classification.done",
    output:
        done = touch("{results}/done/pangenome.mlsa_find_genes.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['mlsa_find_genes.opts'],
        suggestions = "{results}/pangenome_db/gene_classification/mlsa_suggestions.txt",
    benchmark:
        "{results}/benchmarks/pangenome.mlsa_find_genes.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.6
    shell:
        """
        genes=`grep -Ev ':|^#|^$' {params.suggestions} | sed -z 's/\\n/,/g' | sed -e 's/,*$//g' -e 's/_mRNA//g'`
        {pantools} mlsa_find_genes --genes=$genes {params.opts} {params.database}
        """

rule mlsa_concatenate:
    """Step 2/3 of mlsa. Concatenate the gene selection into a single continuous sequence."""
    input:
        "{results}/done/pangenome.mlsa_find_genes.done"
    output:
        done = touch("{results}/done/pangenome.mlsa_concatenate.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['mlsa_concatenate.opts'],
        suggestions = "{results}/pangenome_db/gene_classification/mlsa_suggestions.txt",
    benchmark:
        "{results}/benchmarks/pangenome.mlsa_concatenate.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        genes=`grep -Ev ':|^#|^$' {params.suggestions} | sed -z 's/\\n/,/g' | sed -e 's/,*$//g' -e 's/_mRNA//g'`
        {pantools} mlsa_concatenate -t={threads} -g=$genes {params.opts} {params.database}
        """

rule mlsa:
    """Step 3/3 of mlsa. Run IQ-tree on the concatenated sequences."""
    input:
        "{results}/done/pangenome.mlsa_concatenate.done",
        "{results}/done/pangenome.add_phenotypes.done" if config['phenotypes'] else [],
    output:
        done = touch("{results}/done/pangenome.mlsa.done"),
    params:
        database = "{results}/pangenome_db",
        opts = config['mlsa.opts'],
    benchmark:
        "{results}/benchmarks/pangenome.mlsa.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "{pantools} mlsa -t={threads} {params.opts} {params.database}"