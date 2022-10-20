# Rules to assess synteny between the haplotypes within subsegenomes

rule minimap_toPAF:
    input:
        unpack(get_minimap_input_files)
    output:
        f'{SYNTENY_DIR}/minimap/{{sg}}_hap1_vs_hap2.paf'
    conda: '../envs/synteny.yaml'
    log: LOG_DIR + '/minimap_toPAF/{sg}_minimap_toPAF.log'
    threads: 8
    shell:
        """
        ( minimap2 -cx asm5 {input.hap1} {input.hap2} \
            -t {threads} \
            --cs | sort -k6,6 -k8,8 > {output} ) 2> {log}
        """

rule minimap_toBAM:
    input:
        unpack(get_minimap_input_files)
    output:
        f'{SYNTENY_DIR}/minimap/{{sg}}_hap1_vs_hap2.bam'
    conda: '../envs/synteny.yaml'
    log: LOG_DIR + '/minimap_toBAM/{sg}_minimap_toBAM.log'
    threads: 8
    shell:
        """
        ( minimap2 -ax asm5 {input.hap1} {input.hap2} \
            -t {threads} \
            --cs | 
                samtools sort -@ {threads} \
                    -O bam \
                    -T {wildcards.sg}_sorting.tmp \
                    -o {output} - ) 2> {log}
        """

rule synteny_done:
    input:
        expand(rules.minimap_toPAF.output, sg=['occ', 'pall']),
        expand(rules.minimap_toBAM.output, sg=['occ', 'pall'])
    output:
        f"{SYNTENY_DIR}/synteny.done"
    shell:
        """
        touch {output}
        """


