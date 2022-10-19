# Python functions used through workflow

def get_haplotype_fasta(wildcards):
    path = f"{HAPLOTYPE_FASTA_DIR}/{wildcards.hap}/"
    fasta = glob.glob(path + '*.fasta')[0]
    return fasta
