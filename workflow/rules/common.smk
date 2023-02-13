# Python functions used through workflow

def get_haplotype_fasta(wildcards):
    fasta = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/{wildcards.hap}.fasta")[0]
    return fasta

def get_minimap_hap_vs_hap_input_files(wildcards):
    hap1 = wildcards.hap_comp.split('-')[0]
    hap2 = wildcards.hap_comp.split('-')[1]
    
    if wildcards.ver == 'original':
        hap1_fasta = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/{hap1}.fasta")[0]
        hap2_fasta = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/{hap2}.fasta")[0]
    elif wildcards.ver == 'revised':
        all_fastas = expand(rules.fix_haplotypes.output, hap = HAPS)
        hap1_fasta = [f for f in all_fastas if hap1 in f][0]
        hap2_fasta = [f for f in all_fastas if hap2 in f][0]

    return { 'hap1' : hap1_fasta, 'hap2' : hap2_fasta }
    
def get_star_align_input_files(wildcards):
    star_build = rules.build_star.output
    if wildcards.acc.startswith('TR'):
        R1 = glob.glob(f"{config['kooyers_rnaseq']}/{wildcards.acc}/{wildcards.acc}_*_1.fq.gz")[0]
        R2 = glob.glob(f"{config['kooyers_rnaseq']}/{wildcards.acc}/{wildcards.acc}_*_2.fq.gz")[0]
    else:
        R1 = rules.gzip_fastq.output.R1
        R2 = rules.gzip_fastq.output.R2
    return { 'R1' : R1, 'R2' : R2, 'star_build' : star_build }
        
