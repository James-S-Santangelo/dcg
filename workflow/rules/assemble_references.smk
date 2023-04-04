# Rules to create haploid reference assembly and phased diploid assembly for annotation and publishing.

rule fix_haplotypes:
    input:
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.agp" for hap in HAPS],
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.fasta" for hap in HAPS],
        expand(rules.dotplot_hap_vs_hap.output, ver='original', hap_comp=HAP_VS_HAP)
    output:
        expand(f"{REVISED_HAP_DIR}/{{hap}}_revised.fasta", hap=HAPS)
    log: LOG_DIR + '/notebooks/fix_haplotypes_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/fix_assembly_issues.py.ipynb"

rule TrRvSix_vs_TrRvFive_and_LG:
    input:
        expand(rules.blast_markers.output, hap=HAPS, map_pop=['SG','DG']),
        f"{GENMAP_RESOURCE_DIR}/DG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/DG_genMap.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_genMap.csv",
        TRR_FIVE_SCAFFS,
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS)
    output:
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/BLASThits_linkageGroup_by_subGenome.pdf',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_scaffoldLengths_bySubGenome_and_Hap.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrRv6_longestScaffolds_byLinkageGroup.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeMapping.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeOrientation.pdf',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_chromosomesToReverseComplement.csv',
    log: LOG_DIR + '/notebooks/TrRvSix_vs_TrRvFive_and_LG_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/TrRvSix_vs_TrRvFive_and_LG.r.ipynb"

rule create_reference_assemblies:
    input:
        rules.fix_haplotypes.output,
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeMapping.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_chromosomesToReverseComplement.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_scaffoldLengths_bySubGenome_and_Hap.csv'
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_haploid_reference.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/diploid_reference/TrR_v6_hap1_reference.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/diploid_reference/TrR_v6_hap2_reference.fasta",
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
