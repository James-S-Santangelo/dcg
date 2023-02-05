# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule dotplot_hap_vs_hap:
    input:
        rules.minimap_hap_vs_hap_paf.output
    output:
        f"{QC_DIR}/dotplots/{{ver}}_haps/{{hap_comp}}.pdf"
    conda: '../envs/qc.yaml'
    script:
        "../scripts/r/dotplot_hap_vs_hap.R"

rule qc_done:
    input:
        expand(rules.dotplot_hap_vs_hap.output, ver=['original', 'revised'], hap_comp=HAP_VS_HAP)
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
