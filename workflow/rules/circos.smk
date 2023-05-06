
###############
#### SETUP ####
###############

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
            -w 100000  \
            -s 20000 \
            -i srcwinnum > {output.bed} 2> {log}
        """

rule makeblastdb_fromRef:
    """
    Creates BLAST Database from UTM reference or Griffiths reference (chromosomes only).
    """
    input:
        lambda w: rules.split_fasta_toChroms_andOrganelles.output.chroms if w.ref == 'utm' else rules.TrRvFive_chromsOnly.output.fasta
    output:
        multiext(f'{BLAST_DIR}/blastDBs/{{ref}}/{{ref}}.fasta', '.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto') 
    conda: '../envs/blast.yaml'
    log: LOG_DIR + '/makeblastdb/makeblastdb_{ref}.log'
    params:
        outfile = f'{BLAST_DIR}/blastDBs/{{ref}}/{{ref}}.fasta'
    shell:
        """
        makeblastdb -in {input} \
            -dbtype nucl \
            -parse_seqids \
            -logfile {log} \
            -out {params.outfile}
        """

############################
#### LINKAGE MAP CIRCOS ####
############################

rule blast_linkageMarkers_toRef:
    """
    BLASTs linkage markers from Olsen et al. (2022) F2 mapping populations against UTM or Griggiths references. Used to generate Circos plot with markers connected to both reference 
    """
    input:
        markers = ancient(f'{GENMAP_RESOURCE_DIR}/{{map_pop}}_tags.fa'),
        db = rules.makeblastdb_fromRef.output 
    output:
        f'{BLAST_DIR}/genMap/{{ref}}_{{map_pop}}_marker_blast.txt'
    conda: '../envs/blast.yaml'
    log: LOG_DIR + '/blast_markers/{ref}_{map_pop}_marker_blast.log'
    params:
        outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp'",
        db_base = lambda wildcards, input: os.path.splitext(input.db[0])[0]
    threads: 4
    shell:
        """
        blastn -db {params.db_base} \
            -query {input.markers} \
            -out {output} \
            -outfmt {params.outfmt} \
            -num_threads {threads} \
            -evalue 1e-10 \
            -max_hsps 5 \
            -max_target_seqs 5 2> {log}  
        """

# rule create_karyotype_file_linkageMap:
#     input:
#         rules.get_chrom_lengths.output
#     output:
#         f"{FIGURES_DIR}/circos/data/utm_karyotype.txt"
#     run:
#         if wildcards.ref == 'utm':
#             color = 'blue'
#         else:
#             color = 'orange'
#         with open(output[0], 'w') as fout:
#             with open(input[0], 'r') as fin:
#                 if wildcards.ref == 'utm':
#                     lines = fin.readlines()
#                 else:
#                     lines = reversed(fin.readlines())
#                 for line in lines:
#                     sline = line.strip().split('\t')
#                     fout.write(f"chr\t-\t{sline[0]}\t{sline[0]}\t0\t{sline[1]}\t{color}\n")

##############################
#### UTM_TREP_v1.0 CIRCOS #### 
##############################

#### GC CONTENT ####

rule bedtools_nuc:
    input:
        fasta = rules.split_fasta_toChroms_andOrganelles.output.chroms,
        bed = rules.bedtools_makewindows.output.bed
    output:
        f"{FIGURES_DIR}/circos/gc/utm_windowed_gc_content.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedtools/utm_nuc.log"
    shell:
        """
        bedtools nuc -fi {input.fasta} -bed {input.bed} > {output} 2> {log}
        """

#### GENE AND REPEAT DENSITY ####

rule gff_genesOnly:
    input:
        rules.gff_sort_functional.output 
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/utm_genes.gff"
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
        get_gffToBed_input 
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/utm_{{feat}}.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/utm_{{feat}}_toBed.log"
    shell:
        """
        gff2bed < {input} > {output} 2> {log}
        """

rule chromLengths_toBed:
    input:
        bed = rules.get_chrom_lengths.output
    output:
         f"{PROGRAM_RESOURCE_DIR}/circos/utm_chrom_lengths.bed"
    conda: '../envs/circos.yaml'
    shell:
        """
        awk -vOFS="\t" '{{ print $1, "0", $2; }}' {input.bed} | sort-bed - > {output}
        """

rule bedops:
    input:
        bed = rules.chromLengths_toBed.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/utm_bedops_windows.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/utm_bedops.log"
    shell:
        """
        bedops \
            --chop 100000 \
            --stagger 20000 \
            -x {input.bed} > {output} 2> {log}
        """

rule bedmap_featureCount:
    input:
        win = rules.bedops.output,
        feat_bed = rules.gffToBed.output
    output:
        f"{FIGURES_DIR}/circos/{{feat}}/utm_{{feat}}_count.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/feature_count/utm_{{feat}}_count.log"
    shell:
        """
        bedmap \
            --echo \
            --count \
            --delim '\t' \
            {input.win} {input.feat_bed} > {output} 2> {log}
        """
      
rule circos_done:
    input:
       expand(rules.blast_linkageMarkers_toRef.output, ref = ['utm', 'TrRv5'], map_pop = ['SG', 'DG']) 
    output:
        f"{FIGURES_DIR}/circos/circos.done"
    shell:
        """
        touch {output}
        """
        






