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
        lambda w: rules.index_fasta_chromsOnly.output if w.ref == 'utm' else f"{TRR_FIVE_FASTA}.fai"
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/{{ref}}_chrom_lengths.txt"
    conda: '../envs/circos.yaml'
    shell:
        """
        if [ {wildcards.ref} = 'utm' ]; then
            cut -f1,2 {input} > {output}
        else
            grep '^CM' {input} | cut -f1,2 > {output}
        fi
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
        fasta = lambda w: rules.split_fasta_toChroms_andOrganelles.output.chroms if w.ref == 'utm' else TRR_FIVE_FASTA,
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
        bed = rules.bedtools_makewindows.output
    output:
         f"{PROGRAM_RESOURCE_DIR}/circos/chrom_lengths.bed"
    conda: '../envs/circos.yaml'
    shell:
        """
        awk -vOFS="\t" '{{ print $1, "0", $2; }}' {input.bed} | sort-bed - > {output}
        """

rule gff_genesOnly:
    input:
        rules.gff_sort_functional.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/gene.gff"
    run:
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                lines = fin.readlines()
                for line in lines:
                    if not line.startswith('#'):
                        sline = line.split('\t')
                        if sline[2] == 'gene':
                            fout.write(f"{line}")

rule gffToBed:
    input:
        gff = lambda w: rules.repeat_masker.output.gff if w.feat == 'repeat' else rules.gff_genesOnly.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/{{feat}}.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/{{feat}}_toBed.log"
    shell:
        """
        gff2bed < {input.gff} > {output} 2> {log}
        """

rule bedops:
    input:
        bed = rules.chromLengths_toBed.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/bedops_windows.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/bedops.log"
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
        f"{FIGURES_DIR}/circos/{{feat}}/{{feat}}_count.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/feature_count/{{feat}}_count.log"
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
        ref = lambda w: rules.split_fasta_toChroms_andOrganelles.output.chroms if w.ref == 'utm' else TRR_FIVE_FASTA,
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

rule circos_done:
    input:
        expand(rules.bedtools_nuc.output, ref = ['utm','TrRv5']),
        expand(rules.bedmap_featureCount.output, feat = ['repeat', 'gene']),
        expand(rules.bedtools_multicov.output, ref = ['utm','TrRv5'])
    output:
        f"{FIGURES_DIR}/circos/circos.done"
    shell:
        """
        touch {output}
        """
        






