# Rules to create haploid reference assembly and phased diploid assembly for annotation and publishing.

rule linkage_analysis:
    input:
        expand(rules.blast_markers.output, hap=HAPS, map_pop=['SG','DG']),
        f"{GENMAP_RESOURCE_DIR}/DG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/DG_genMap.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_genMap.csv",
        config['TrR_v5_scaffs'],
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

rule reference_assemblies_done:
    input:
        rules.linkage_analysis.output
    output:
        f'{REFERENCE_ASSEMBLIES_DIR}/reference_assemblies.done'
    shell:
        """
        echo 'Done'
        """
