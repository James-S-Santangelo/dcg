# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule dotplot_hap_vs_hap:
    input:
        rules.minimap_hap_vs_hap_paf.output
    output:
        f"{QC_DIR}/dotplots/{{ver}}_haps/{{hap_comp}}.pdf"
    conda: '../envs/qc.yaml'
    script:
        "../scripts/r/dotplot_hap_vs_hap.R"

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
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrRv5_and_TrRv6_longestScaffolds_byLinkageGroup.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeMapping.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeOrientation.pdf',
    log: LOG_DIR + '/notebooks/TrRvSix_vs_TrRvFive_and_LG_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/TrRvSix_vs_TrRvFive_and_LG.r.ipynb"

rule qc_done:
    input:
        expand(rules.dotplot_hap_vs_hap.output, ver=['original', 'revised'], hap_comp=HAP_VS_HAP),
        rules.TrRvSix_vs_TrRvFive_and_LG.output
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
