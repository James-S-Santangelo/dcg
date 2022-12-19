# Python functions used through workflow

def get_haplotype_fasta(wildcards):
    path = f"{HAPLOTYPE_FASTA_DIR}/{wildcards.hap}/"
    if wildcards.hap == 'pall1':
        fasta = glob.glob(path + '*.fasta')[0]
    elif wildcards.hap == 'pall2':
        fasta = glob.glob(path + 'new_final_assembly.fasta')[0]
    else:
        fasta = glob.glob(path + 'review_final_assembly.fasta')[0]
    return fasta

def get_minimap_input_files(wildcards):
    if wildcards.sg == 'occ':
        hap1 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/occ1/review_final_assembly.fasta")[0]
        hap2 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/occ2/review_final_assembly.fasta")[0]
    else:
        hap1 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/pall1/*.fasta")[0]
        hap2 = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/pall2/new_final_assembly.fasta")[0]
    return { 'hap1' : hap1, 'hap2' : hap2 }
    
def get_star_align_input_files(wildcards):
    star_build = rules.build_star.output
    if wildcards.acc.startswith('TR'):
        R1 = glob.glob(f"{config['kooyers_rnaseq']}/{wildcards.acc}/{wildcards.acc}_*_1.fq.gz")[0]
        R2 = glob.glob(f"{config['kooyers_rnaseq']}/{wildcards.acc}/{wildcards.acc}_*_2.fq.gz")[0]
    else:
        R1 = rules.fasterq_dump.output.R1
        R2 = rules.fasterq_dump.output.R2
    return { 'R1' : R1, 'R2' : R2, 'star_build' : star_build }
        
