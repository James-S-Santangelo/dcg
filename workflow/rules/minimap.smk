# Rules to assess minimap between the haplotypes within subsegenomes

rule minimap_hap_vs_hap_paf:
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

rule minimap_done:
    input:
        expand(rules.minimap_hap_vs_hap_paf.output, ver=['original', 'revised'], hap_comp=HAP_VS_HAP),
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS)
    output:
        f"{MINIMAP_DIR}/minimap.done"
    shell:
        """
        touch {output}
        """


