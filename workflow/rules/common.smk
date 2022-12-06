# Python functions used through workflow

def get_haplotype_fasta(wildcards):
    path = f"{HAPLOTYPE_FASTA_DIR}/{wildcards.hap}/"
    if wildcards.hap == 'pall1':
        fasta = glob.glob(path + '*.fasta')[0]
    else:
        fasta = glob.glob(path + 'new_final_assembly.fasta')[0]
    return fasta

def get_minimap_input_files(wildcards):
    if wildcards.sg == 'occ':
        hap1 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/occ1/new_final_assembly.fasta")[0]
        hap2 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/occ2/new_final_assembly.fasta")[0]
    else:
        hap1 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/pall1/*.fasta")[0]
        hap2 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/pall2/new_final_assembly.fasta")[0]
    return { 'hap1' : hap1, 'hap2' : hap2 }
        
