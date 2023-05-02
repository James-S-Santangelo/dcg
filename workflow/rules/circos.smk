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

rule bedtools_makewindows:
    input:
        fai = rules.index_fasta_chromsOnly.output
    output:
        bed = f"{FIGURES_DIR}/circos/gc/gc_content_windows.bed",
        chr_len = f"{PROGRAM_RESOURCE_DIR}/circos/chrom_lengths.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedtools/makewindows.log"
    shell:
        """
        cut -f1,2 {input.fai} > {output.chr_len}
        bedtools makewindows \
            -g {output.chr_len} \
            -w 1000000  \
            -s 200000 \
            -i srcwinnum > {output.bed} 2> {log}
        """

rule bedtools_nuc:
    input:
        fasta = rules.split_fasta_toChroms_andOrganelles.output.chroms,
        bed = rules.bedtools_makewindows.output.bed
    output:
        f"{FIGURES_DIR}/circos/gc/windowed_gc_content.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedtools/nuc.log"
    shell:
        """
        bedtools nuc -fi {input.fasta} -bed {input.bed} > {output} 2> {log}
        """

rule chromLengths_toBed:
    input:
        bed = rules.bedtools_makewindows.output.chr_len
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
        
