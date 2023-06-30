"""
Validate paths in the provided location files.
If paths refer to the root directory of a downloaded set within the resources directory, 
they will be adjusted to go from the root of the pipeline.
"""

rule validate_genome_paths:
    input:
        config['genomes']
    output:
        done = touch("{results}/done/validate.genomes.done")
    script:
        "../scripts/validate_locations.py"

rule validate_proteome_paths:
    input:
        config['proteomes']
    output:
        done = touch("{results}/done/validate.proteomes.done")
    script:
        "../scripts/validate_locations.py"

rule validate_annotation_paths:
    input:
        config['annotations']
    output:
        done = touch("{results}/done/validate.annotations.done")
    script:
        "../scripts/validate_locations.py"

rule validate_function_paths:
    input:
        config['functions']
    output:
        done = touch("{results}/done/validate.functions.done")
    script:
        "../scripts/validate_locations.py"

rule validate_pav_paths:
    input:
        config['pav']
    output:
        done = touch("{results}/done/validate.pavs.done")
    script:
        "../scripts/validate_locations.py"

rule validate_vcf_paths:
    input:
        config['vcf']
    output:
        done = touch("{results}/done/validate.variants.done")
    script:
        "../scripts/validate_locations.py"