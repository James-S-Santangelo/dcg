# Rules to assess minimap between the haplotypes within subsegenomes

rule minimap_hap_vs_hap_paf:
    """
    Maps each haplotype against the other using Minimap. This is done for the raw original haplotypes, and the revised haplotypes following manual fixes. 

    The alignments for the raw haplotypes were used for QC, identifying misallignments, and identifying scaffolds to fill gaps in the original haplotypes. 
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

rule combine_haps:
    """
    Rule to separately combine Pall and Occ Hap1, and Pall and Occ Hap2 for global Minimap alignment
    """
    input:
        occ_hap1 = lambda w: f"{HAPLOTYPE_FASTA_DIR}/occ1.fasta" if w.ver == 'original' else f"{REVISED_HAP_DIR}/occ1_revised.fasta",
        occ_hap2 = lambda w: f"{HAPLOTYPE_FASTA_DIR}/occ2.fasta" if w.ver == 'original' else f"{REVISED_HAP_DIR}/occ2_revised.fasta",
        pall_hap1 = lambda w: f"{HAPLOTYPE_FASTA_DIR}/pall1.fasta" if w.ver == 'original' else f"{REVISED_HAP_DIR}/pall1_revised.fasta",
        pall_hap2 = lambda w: f"{HAPLOTYPE_FASTA_DIR}/pall2.fasta" if w.ver == 'original' else f"{REVISED_HAP_DIR}/pall2_revised.fasta"
    output:
        hap1 = f"{PROGRAM_RESOURCE_DIR}/minimap/{{ver}}/hap1.fasta",
        hap2 = f"{PROGRAM_RESOURCE_DIR}/minimap/{{ver}}/hap2.fasta"
    conda: '../envs/minimap.yaml'
    shell:
        """
        sed 's/__.*//' {input.occ_hap1} | sed 's/>Scaffold/>o1S/' > {wildcards.ver}tmp_occ1.fasta
        sed 's/__.*//' {input.occ_hap2} | sed 's/>Scaffold/>o2S/' > {wildcards.ver}tmp_occ2.fasta
        sed 's/__.*//' {input.pall_hap1} | sed 's/>Scaffold/>p1S/' > {wildcards.ver}tmp_pall1.fasta
        sed 's/__.*//' {input.pall_hap2} | sed 's/>Scaffold/>p2S/' > {wildcards.ver}tmp_pall2.fasta

        awk 1 {wildcards.ver}tmp_occ1.fasta {wildcards.ver}tmp_pall1.fasta > {wildcards.ver}tmp1.fasta
        awk 1 {wildcards.ver}tmp_occ2.fasta {wildcards.ver}tmp_pall2.fasta > {wildcards.ver}tmp2.fasta

        samtools faidx {wildcards.ver}tmp1.fasta
        samtools faidx {wildcards.ver}tmp2.fasta

        cat {wildcards.ver}tmp1.fasta.fai | awk '$2 > 3000000 {{print $1}}' > {wildcards.ver}tmp1_scaffs.txt
        cat {wildcards.ver}tmp2.fasta.fai | awk '$2 > 3000000 {{print $1}}' > {wildcards.ver}tmp2_scaffs.txt

        samtools faidx {wildcards.ver}tmp1.fasta -r {wildcards.ver}tmp1_scaffs.txt > {output.hap1}
        samtools faidx {wildcards.ver}tmp2.fasta -r {wildcards.ver}tmp2_scaffs.txt > {output.hap2}

        rm {wildcards.ver}tmp*
        """

rule minimap_hapsCombined_paf:
    """
    Maps Hap1 against Hap2, rather than pairwise. This is done for the raw original haplotypes, and the revised haplotypes following manual fixes. 

    The alignments for the raw haplotypes were used for QC, identifying misallignments, and identifying scaffolds to fill gaps in the original haplotypes. Used to generate dotplots. 
    """
    input:
        hap1 = rules.combine_haps.output.hap1,
        hap2 = rules.combine_haps.output.hap2
    output:
        f'{MINIMAP_DIR}/haps_combined/{{ver}}/allHap1_vs_allHap2.paf'
    conda: '../envs/minimap.yaml'
    log: LOG_DIR + '/minimap/hapsCombined_{ver}_paf.log'
    threads: 16
    shell:
        """
        ( minimap2 -cx asm5 {input.hap2} {input.hap1} \
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
    Generate simple dotplot of alignments of pairwise haplotype alignments.
    """
    input:
        rules.minimap_hap_vs_hap_paf.output
    output:
        f"{FIGURES_DIR}/dotplots/{{ver}}_haps/{{hap_comp}}.pdf"
    conda: '../envs/qc.yaml'
    script:
        "../scripts/r/dotplot_hap_vs_hap.R"

rule dotplot_haps_combined:
    """
    Generate simple dotplot of alignments of one haplotype mapped against the other.
    """
    input:
        rules.minimap_hapsCombined_paf.output
    output:
        f"{FIGURES_DIR}/dotplots/haps_combined/{{ver}}_haps.pdf"
    conda: '../envs/qc.yaml'
    script:
        "../scripts/r/dotplot_hap_vs_hap.R"

rule minimap_done:
    input:
        expand(rules.dotplot_hap_vs_hap.output, ver=['original', 'revised'], hap_comp=HAP_VS_HAP),
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS),
        expand(rules.dotplot_haps_combined.output, ver = ['original', 'revised'])
    output:
        f"{MINIMAP_DIR}/minimap.done"
    shell:
        """
        touch {output}
        """


