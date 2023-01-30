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

rule create_reference_assemblies:
    input:
       f'{REFERENCE_ASSEMBLIES_DIR}/resources/TrR_v5_to_v6_chromosomeMapping.csv',
        expand(rules.minimap_hap_vs_hap_paf.output, hap_comp=HAP_VS_HAP)
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_haploid_reference.fasta"
    conda:
        '../envs/notebooks.yaml'
    params:
        HAPLOTYPE_FASTA_DIR
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
