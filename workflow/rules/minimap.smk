# Rules to assess minimap between the haplotypes within subsegenomes

rule minimap_hap_vs_hap_paf:
    """
    Maps each haplotype against the other using Minimap. This is done for the raw original haplotypes, and the revised haplotypes following manual fixes. 

    The alignments for the raw haplotypes were used for QC, identifying misallignments, and identifying scaffolds to fill gaps in the original haplotypes. Used to generate dotplots in figure S1. 
    """
    input:
        unpack(get_minimap_hap_vs_hap_input_files)
    output:
        f'{MINIMAP_DIR}/hap_vs_hap/{{ver}}/{{hap_comp}}.paf'
    conda: '../envs/minimap.yaml'
    log: LOG_DIR + '/minimap/{hap_comp}_{ver}_paf.log'
    threads: 8
    shell:
        """
        ( minimap2 -cx asm5 {input.hap1} {input.hap2} \
            -t {threads} \
            --cs | sort -k6,6 -k8,8 > {output} ) 2> {log}
        """

rule minimap_hap_vs_TrRvFive:
    """
    Maps each haplotype against the previous white clover assembly using Minimap. Used for assessing correspondance between previous chromosomes and chromosomes in the new haplotype assemblies. 
    """
    input:
        hap = get_haplotype_fasta,
        TrR_five = TRR_FIVE_FASTA 
    output:
        f"{MINIMAP_DIR}/hap_vs_TrRv5/{{hap}}_vs_TrRv5.paf"
    conda: '../envs/minimap.yaml'
    log: LOG_DIR + '/minimap/{hap}_vs_TrRv5.log'
    threads: 8
    shell:
        """
        ( minimap2 -cx asm5 {input.hap} {input.TrR_five} \
            -t {threads} \
            --cs | sort -k6,6 -k8,8 > {output} ) 2> {log}
        """

rule dotplot_hap_vs_hap:
    """
    Generate simple dotplot of alignments of one haplotype mapped against the other.
    """
    input:
        rules.minimap_hap_vs_hap_paf.output
    output:
        f"{FIGURES_DIR}/dotplots/{{ver}}_haps/{{hap_comp}}.pdf"
    conda: '../envs/qc.yaml'
    script:
        "../scripts/r/dotplot_hap_vs_hap.R"

rule minimap_done:
    input:
        expand(rules.dotplot_hap_vs_hap.output, ver=['original', 'revised'], hap_comp=HAP_VS_HAP),
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS)
    output:
        f"{MINIMAP_DIR}/minimap.done"
    shell:
        """
        touch {output}
        """


