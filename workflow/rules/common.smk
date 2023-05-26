# Python functions used through workflow

def get_haplotype_fasta(wildcards):
    """
    Simple Python Function that retrieves the raw FASTA file for a given haplotype
    """
    fasta = glob.glob(f"{HAPLOTYPE_FASTA_DIR}/{wildcards.hap}.fasta")[0]
    return fasta


def get_minimap_hap_vs_hap_input_files(wildcards):
    """
    Retrieves input files for rule mapping each haplotype against the other using Minimap. 
    """
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


def get_create_subgenome_fasta_input(wildcards):
    all_fasta = expand(rules.split_chromosomal_fasta.output, chrom=CHROMOSOMES)
    sg_fastas = [fa for fa in all_fasta if wildcards.sg in os.path.basename(fa)]
    return sg_fastas
