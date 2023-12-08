"""
Contains setup logic and common functions for other Snakemake rules.
"""

import tempfile

# Set pantools command with Java heap space.
if config['jar']:
    pantools = "java -Xms{}g -Xmx{}g -jar {}".format(config['Xms'], config['Xmx'], config['jar'])
else:
    pantools = "pantools -Xms{}g -Xmx{}g".format(config['Xms'], config['Xmx'])

if config['nice']: pantools = f"nice -n {config['nice']} {pantools}"

# Set the temporary directory based on system default or given path.
temp = config['scratch'] if config['scratch'] else tempfile.gettempdir()
temp_dir = tempfile.TemporaryDirectory(dir = temp).name

# Gives needed input for functions that use proteins depending on db type.
def proteins_done(db_type):
    if db_type == 'pangenome':
        return "{results}/done/pangenome.add_annotations.done"
    else: 
        return "{results}/done/panproteome.build_panproteome.done"

def msa_done(db_type):
    if config['vcf']:
        return "{results}/done/{type}.msa_variants.done"
    if db_type == 'pangenome':
        return "{results}/done/{type}.msa_nucleotide.done"
    else:
        return "{results}/done/{type}.msa_protein.done"

def construction_done(db_type):
    done_files = []
    if db_type == 'pangenome':
        done_files.append("{results}/done/{type}.build_pangenome.done")
        done_files.append("{results}/done/{type}.add_annotations.done")
        if config['vcf']: done_files.append("{results}/done/{type}.add_variants.done")
        if config['phasing']: done_files.append("{results}/done/{type}.add_phasing.done")
    else:
        done_files.append("{results}/done/{type}.build_panproteome.done")
    if config['functions']: done_files.append("{results}/done/{type}.add_functions.done")
    if config['pav']: done_files.append("{results}/done/{type}.add_pavs.done")
    if config['phenotypes']: done_files.append("{results}/done/{type}.add_phenotypes.done")
    return done_files