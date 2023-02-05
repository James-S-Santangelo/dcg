# Rules to create haploid reference assembly and phased diploid assembly for annotation and publishing.

rule linkage_analysis:
    input:
        expand(rules.blast_markers.output, hap=HAPS, map_pop=['SG','DG']),
        f"{GENMAP_RESOURCE_DIR}/DG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/DG_genMap.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_genMap.csv",
        TRR_FIVE_SCAFFS,
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS)
    output:
        f'{REFERENCE_ASSEMBLIES_DIR}/resources/BLASThits_LGxSG.pdf',
        f'{REFERENCE_ASSEMBLIES_DIR}/resources/scaffoldLengths_bySG_and_Hap_TrR_v6.csv',
        f'{REFERENCE_ASSEMBLIES_DIR}/resources/longestScaffolds_byLG_TrR_v5_and_v6.csv',
        f'{REFERENCE_ASSEMBLIES_DIR}/resources/TrR_v5_to_v6_chromosomeMapping.csv',
        f'{REFERENCE_ASSEMBLIES_DIR}/resources/TrR_v5_to_v6_chromosomeOrientation.pdf',
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/linkage_analysis.r.ipynb"

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
        f'{REFERENCE_ASSEMBLIES_DIR}/resources/TrR_v5_to_v6_chromosomeMapping.csv',
        expand(rules.minimap_hap_vs_hap_paf.output, ver='original', hap_comp=HAP_VS_HAP)
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_haploid_reference.fasta"
    conda:
        '../envs/notebooks.yaml'
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
