# Rules to create haploid reference assembly and phased diploid assembly for annotation and publishing.

rule fix_haplotype_issues:
    input:
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.agp" for hap in HAPS],
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.fasta" for hap in HAPS],
        expand(rules.dotplot_hap_vs_hap.output, ver='original', hap_comp=HAP_VS_HAP)
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/revised_haplotypes/occ1_revised.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/revised_haplotypes/occ2_revised.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/revised_haplotypes/pall1_revised.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/revised_haplotypes/pall2_revised.fasta"
    log: LOG_DIR + '/notebooks/fix_haplotype_issues_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/fix_assembly_issues.py.ipynb"

rule create_reference_assemblies:
    input:
        rules.fix_haplotype_issues.output,
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeMapping.csv',
        expand(rules.minimap_hap_vs_hap_paf.output, ver='original', hap_comp=HAP_VS_HAP)
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_haploid_reference.fasta"
    log: LOG_DIR + '/notebooks/generate_assemblies_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/generate_assemblies.py.ipynb"

rule reference_assemblies_done:
    input:
        rules.create_reference_assemblies.output
    output:
        f'{REFERENCE_ASSEMBLIES_DIR}/reference_assemblies.done'
    shell:
        """
        touch {output}
        """
