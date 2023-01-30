# Rules to assess synteny between the haplotypes within subsegenomes

rule minimap_hap_vs_hap_paf:
    input:
        unpack(get_minimap_hap_vs_hap_input_files)
    output:
        f'{SYNTENY_DIR}/minimap/hap_vs_hap/paf/{{hap_comp}}.paf'
    conda: '../envs/synteny.yaml'
    log: LOG_DIR + '/minimap/minimap_{hap_comp}_paf.log'
    threads: 8
    shell:
        """
        ( minimap2 -cx asm5 {input.hap1} {input.hap2} \
            -t {threads} \
            --cs | sort -k6,6 -k8,8 > {output} ) 2> {log}
        """


rule minimap_hap_vs_hap_bam:
    input:
        unpack(get_minimap_hap_vs_hap_input_files)
    output:
        f'{SYNTENY_DIR}/minimap/hap_vs_hap/bam/{{hap_comp}}.bam'
    conda: '../envs/synteny.yaml'
    log: LOG_DIR + '/minimap/minimap_{hap_comp}_log.log'
    threads: 8
    shell:
        """
        ( minimap2 -ax asm5 {input.hap1} {input.hap2} \
            -t {threads} | samtools sort -@ {threads} - -o {output} &&\
            samtools index {output} ) 2> {log}
        """

rule minimap_hap_vs_TrRvFive:
    input:
        hap = get_haplotype_fasta,
        TrR_five = TRR_FIVE_FASTA 
    output:
        f"{SYNTENY_DIR}/minimap/hap_vs_TrRv5/{{hap}}_vs_TrRv5.paf"
    conda: '../envs/synteny.yaml'
    log: LOG_DIR + '/minimap/{hap}_vs_TrRv5.log'
    threads: 8
    shell:
        """
        ( minimap2 -cx asm5 {input.hap} {input.TrR_five} \
            -t {threads} \
            --cs | sort -k6,6 -k8,8 > {output} ) 2> {log}
        """
        

rule synteny_done:
    input:
        expand(rules.minimap_hap_vs_hap_paf.output, hap_comp=HAP_VS_HAP),
        expand(rules.minimap_hap_vs_hap_bam.output, hap_comp=HAP_VS_HAP),
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS)
    output:
        f"{SYNTENY_DIR}/synteny.done"
    shell:
        """
        touch {output}
        """


