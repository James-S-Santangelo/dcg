# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies


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
