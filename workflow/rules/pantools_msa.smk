"""
Run PanTools subcommands for pangenome multiple sequence alignment.
Contains the following subcommands:

msa
"""

rule msa_variants:
    """Create multiple sequence alignments including variants."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_variants.done",
        "{results}/done/{type}.add_functions.done" if config['functions'] else [],
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/done/{type}.add_pavs.done" if config['pav'] else [],
        "{results}/{type}_db/homology_selection.txt" if config['gene_selection'] else []
    output:
        done = touch("{results}/done/{type}.msa_variants.done")
    params:
        database = "{results}/{type}_db",
        opts = config['msa.opts'],
        homology = "-H={results}/{type}_db/homology_selection.txt" if config['gene_selection'] else "",
        pavs = "--pavs" if config['pav'] else ""
    benchmark:
        "{results}/benchmarks/{type}.msa.txt"
    conda:
        "../envs/pantools.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        """
        {pantools} msa \
            --align-variants \
            {params.pavs} \
            -t={threads} \
            {params.opts} \
            {params.homology} \
            {params.database}"""

rule msa_proteins:
    """Create multiple sequence alignments on protein sequences."""
    input:
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_functions.done" if config['functions'] else [],
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/{type}_db/homology_selection.txt" if config['gene_selection'] else []
    output:
        done = touch("{results}/done/{type}.msa_protein.done")
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
        """
        {pantools} msa \
            --align-protein \
            --no-trimming \
            -t={threads} \
            {params.opts} \
            {params.homology} \
            {params.database}
        """


rule msa_nucleotides:
    """Create multiple sequence alignments on nucleotide sequences."""
    input:
        "{results}/done/{type}.msa_protein.done",
        "{results}/done/{type}.grouping.done",
        "{results}/done/{type}.add_functions.done" if config['functions'] else [],
        "{results}/done/{type}.add_phenotypes.done" if config['phenotypes'] else [],
        "{results}/{type}_db/homology_selection.txt" if config['gene_selection'] else []
    output:
        done = touch("{results}/done/{type}.msa_nucleotide.done")
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
        """
        {pantools} msa \
            --align-nucleotide \
            --trim-using-proteins \
            -t={threads} \
            {params.opts} \
            {params.homology} \
            {params.database}
        """
