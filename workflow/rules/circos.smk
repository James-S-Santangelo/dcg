rule TrRvFive_chromsOnly:
    input:
       TRR_FIVE_FASTA
    output:
        fasta = f"{FIGURES_DIR}/circos/TrRv5/TrRv5_chroms.fasta",
        fai = f"{FIGURES_DIR}/circos/TrRv5/TrRv5_chroms.fasta.fai"
    conda: '../envs/circos.yaml'
    params:
        chroms = "CM019101.1 CM019102.1 CM019103.1 CM019104.1 CM019105.1 CM019106.1 CM019107.1 CM019108.1 CM019109.1 CM019110.1 CM019111.1 CM019112.1 CM019113.1 CM019114.1 CM019115.1 CM019116.1"
    shell:
        """
        samtools faidx {input} {params.chroms} > {output.fasta}
        samtools faidx {output.fasta}
        bwa index {output.fasta}
        """

rule index_fasta_chromsOnly:
    input:
        rules.split_fasta_toChroms_andOrganelles.output.chroms
    output:
       f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/split/{ASSEMBLY_NAME}_chromsOnly.fasta.fai"
    conda: '../envs/circos.yaml'
    shell:
        """
        samtools faidx {input}
        """

rule get_chrom_lengths:
    input:
        lambda w: rules.index_fasta_chromsOnly.output if w.ref == 'utm' else rules.TrRvFive_chromsOnly.output.fai
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/{{ref}}_chrom_lengths.txt"
    conda: '../envs/circos.yaml'
    shell:
        """
        cut -f1,2 {input} > {output}
        """

rule bedtools_makewindows:
    input:
        chr_len = rules.get_chrom_lengths.output
    output:
        bed = f"{FIGURES_DIR}/circos/{{ref}}_windows.bed",
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedtools/{{ref}}_makewindows.log"
    shell:
        """
        bedtools makewindows \
            -g {input.chr_len} \
            -w 1000000  \
            -s 200000 \
            -i srcwinnum > {output.bed} 2> {log}
        """

rule bedtools_nuc:
    input:
        fasta = lambda w: rules.split_fasta_toChroms_andOrganelles.output.chroms if w.ref == 'utm' else rules.TrRvFive_chromsOnly.output.fasta,
        bed = rules.bedtools_makewindows.output.bed
    output:
        f"{FIGURES_DIR}/circos/gc/{{ref}}_windowed_gc_content.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedtools/{{ref}}_nuc.log"
    shell:
        """
        bedtools nuc -fi {input.fasta} -bed {input.bed} > {output} 2> {log}
        """

rule chromLengths_toBed:
    input:
        bed = rules.get_chrom_lengths.output
    output:
         f"{PROGRAM_RESOURCE_DIR}/circos/{{ref}}_chrom_lengths.bed"
    conda: '../envs/circos.yaml'
    shell:
        """
        awk -vOFS="\t" '{{ print $1, "0", $2; }}' {input.bed} | sort-bed - > {output}
        """

rule gff_transcriptsOnly:
    input:
        lambda w: rules.gff_sort_functional.output if w.ref == 'utm' else TRR_FIVE_GTF
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/{{ref}}_transcripts.gff"
    run:
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                lines = fin.readlines()
                for line in lines:
                    if not line.startswith('#'):
                        sline = line.split('\t')
                        if wildcards.ref == 'utm':
                            if sline[2] == 'mRNA':
                                fout.write(f"{line}")
                        elif wildcards.ref == 'TrRv5':
                            if sline[2] == 'transcript':
                                fout.write(f"{line}")

rule TrRvFive_repeatMask:
    """
    Softmask repeats in the white clover haploid mapping reference assembly.
    """
    input:
        lib = rules.merge_repeat_databases.output,
        fasta = rules.TrRvFive_chromsOnly.output.fasta
    output:
        fasta = f"{PROGRAM_RESOURCE_DIR}/TrR5_repeat_mask/TrRv5_softMasked.fasta",
        cat = f"{PROGRAM_RESOURCE_DIR}/TrR5_repeat_mask/TrRv5_repeatMasker.cat.gz",
        out = f"{PROGRAM_RESOURCE_DIR}/TrR5_repeat_mask/TrRv5_repeatMasker.out",
        gff = f"{PROGRAM_RESOURCE_DIR}/TrR5_repeat_mask/TrRv5_repeatMasker.gff",
        stats = f"{PROGRAM_RESOURCE_DIR}/TrR5_repeat_mask/TrRv5_repeatMasker.tbl"
    threads: 36
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/TrRvFive_repeat/TrRvFive_repeat.log'
    params:
        outdir = f"{PROGRAM_RESOURCE_DIR}/TrR5_repeat_mask" 
    shell:
        """
        ( RepeatMasker -e ncbi \
            -pa {threads} \
            -nolow \
            -xsmall \
            -gc 33 \
            -lib {input.lib} \
            -dir {params.outdir} \
            -gff {input.fasta} &&
            
            mv {params.outdir}/*.fasta.masked {output.fasta}
            mv {params.outdir}/*.cat.gz {output.cat}
            mv {params.outdir}/*.out {output.out}
            mv {params.outdir}/*.gff {output.gff}
            mv {params.outdir}/*.tbl {output.stats} ) &> {log}
        """

rule gffToBed:
    input:
        get_gffToBed_input 
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/{{ref}}_{{feat}}.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/{{ref}}_{{feat}}_toBed.log"
    shell:
        """
        gff2bed < {input} > {output} 2> {log}
        """

rule bedops:
    input:
        bed = rules.chromLengths_toBed.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/{{ref}}_bedops_windows.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/{{ref}}_bedops.log"
    shell:
        """
        bedops \
            --chop 1000000 \
            --stagger 200000 \
            -x {input.bed} > {output} 2> {log}
        """

rule bedmap_featureCount:
    input:
        win = rules.bedops.output,
        feat_bed = rules.gffToBed.output
    output:
        f"{FIGURES_DIR}/circos/{{feat}}/{{ref}}_{{feat}}_count.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/feature_count/{{ref}}_{{feat}}_count.log"
    shell:
        """
        bedmap \
            --echo \
            --count \
            --delim '\t' \
            {input.win} {input.feat_bed} > {output} 2> {log}
        """
      
rule bwa_index_ref:
    """
    Index reference with bwa to get ready for read mapping
    """
    input:
        rules.split_fasta_toChroms_andOrganelles.output.chroms
    output:
        multiext(f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/split/{ASSEMBLY_NAME}_chromsOnly.fasta", '.amb', '.ann', '.bwt', '.pac', '.sa')
    conda: '../envs/circos.yaml'
    log: f'{LOG_DIR}/bwa/bwa_index_ref.log'
    shell:
        """
        bwa index {input} 2> {log}
        """

rule bwa_map_paired:
    """
    Map trimmed paired reads to the reference genome using bwa mem. Output as BAM
    """
    input:
        r1 = READ1,
        r2 = READ2,
        ref = lambda w: rules.split_fasta_toChroms_andOrganelles.output.chroms if w.ref == 'utm' else rules.TrRvFive_chromsOnly.output.fasta, 
        idx = rules.bwa_index_ref.output
    output:
        f"{FIGURES_DIR}/circos/mapping/{{ref}}.bam" 
    conda: '../envs/circos.yaml'
    log: f'{LOG_DIR}/bwa/{{ref}}_bwa_mem.log'
    threads: 12
    shell:
        """
        ( bwa mem -t {threads} {input.ref} {input.r1} {input.r2} {params} |\
            samtools view -hb -o {output} - ) 2> {log}
        """ 

rule sort_and_index_bam:
    input:
        rules.bwa_map_paired.output
    output:
        bam = f"{FIGURES_DIR}/circos/mapping/{{ref}}_sorted.bam",
        idx = f"{FIGURES_DIR}/circos/mapping/{{ref}}_sorted.bam.bai"
    conda: '../envs/circos.yaml'
    threads: 6
    shell:
        """
        samtools sort -@ {threads} {input} > {output.bam}
        samtools index {output.bam}
        """

rule bedtools_multicov:
    input:
        idx = rules.sort_and_index_bam.output.idx,
        bed = rules.bedtools_makewindows.output.bed,
        bam = rules.sort_and_index_bam.output.bam
    output:
        f"{FIGURES_DIR}/circos/mapping/{{ref}}_windowedCoverage.txt"
    conda: '../envs/circos.yaml'
    log: f'{LOG_DIR}/bedtools/{{ref}}_multicov.log'
    shell:
        """
        bedtools multicov -bed {input.bed} -bams {input.bam} > {output} 2> {log}
        """

rule windowed_MQ:
    input:
        bed = rules.bedtools_makewindows.output,
        bam = rules.sort_and_index_bam.output.bam
    output:
        f"{FIGURES_DIR}/circos/mapping/{{ref}}_windowedMQ.txt"
    conda: '../envs/circos.yaml'
    log: f'{LOG_DIR}/bedtools/{{ref}}_windowedMQ.log'
    shell:
        """
        bedtools map -a {input.bed} \
            -b <(bedtools bamtobed -i {input.bam}) \
            -c 5 -o mean > {output} 2> {log}
        """

rule minimap_utm_vs_TrRvFive:
    """
    Maps the new assembly against the old Griffiths assembly to generate alignments for Circos 
    """
    input:
        utm = rules.split_fasta_toChroms_andOrganelles.output.chroms,
        TrRvFive = rules.TrRvFive_chromsOnly.output.fasta
    output:
        f"{FIGURES_DIR}/circos/minimap/utm_vs_TrRv5.paf"
    conda: '../envs/minimap.yaml'
    log: LOG_DIR + '/minimap/utm_vs_TrRv5.log'
    threads: 8
    shell:
        """
        ( minimap2 -cx asm5 {input.TrRvFive} {input.utm} \
            -t {threads} \
            --cs | sort -k6,6 -k8,8 > {output} ) 2> {log}
        """

rule cutN:
    input:
        ref = lambda w: rules.split_fasta_toChroms_andOrganelles.output.chroms if w.ref == 'utm' else rules.TrRvFive_chromsOnly.output.fasta, 
    output:
        f"{FIGURES_DIR}/circos/cutN/{{ref}}_Ns.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/cutN/{{ref}}_cutN.log"
    shell:
        """
        seqtk cutN -gp10000000 -n1 {input.ref} > {output} 2> {log}
        """

rule circos_done:
    input:
        expand(rules.bedtools_nuc.output, ref = ['utm','TrRv5']),
        expand(rules.bedmap_featureCount.output, feat = ['repeat', 'transcripts'], ref = ['utm','TrRv5']),
        expand(rules.bedtools_multicov.output, ref = ['utm','TrRv5']),
        expand(rules.windowed_MQ.output, ref = ['utm','TrRv5']),
        rules.minimap_utm_vs_TrRvFive.output,
        expand(rules.cutN.output, ref = ['utm','TrRv5'])
    output:
        f"{FIGURES_DIR}/circos/circos.done"
    shell:
        """
        touch {output}
        """
        






