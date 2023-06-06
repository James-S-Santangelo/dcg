
###############
#### SETUP ####
###############

rule TrRvFive_chromsOnly:
    """
    Subset previosu Griffiths assembly into new fasta containing only the chromosomes
    """
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
    """
    Index chromosome-only Griffiths fasta
    """
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
    """
    Get chromosome/scaffold lengths for either Griffiths or new UTM assembly
    """
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
    """
    Use bedtools to create 500 Kb windows with a 100 Kb step
    """
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
            -w 500000  \
            -s 100000 \
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

rule create_karyotype_file_ref:
    """
    Create Circos-formatted karyotype files for Griffiths and UTM reference
    """
    input:
        rules.get_chrom_lengths.output
    output:
        f"{FIGURES_DIR}/circos/{{ref}}_karyotype.txt"
    run:
        if wildcards.ref == 'utm':
            color = "mygreen"
        else:
            color = "myblue"
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                if wildcards.ref == 'utm':
                    lines = reversed(fin.readlines())
                else:
                    lines = reversed(fin.readlines())
                for line in lines:
                    sline = line.strip().split('\t')
                    fout.write(f"chr\t-\t{sline[0]}\t{sline[0]}\t0\t{sline[1]}\t{color}\n")

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

rule generate_circos_genMap_links:
    """
    Create Circos-formatted link and highlight files for linkage markers and thei physical positions in both reference genomes
    """
    input:
        expand(rules.blast_linkageMarkers_toRef.output, map_pop=['SG', 'DG'], ref=['utm', 'TrRv5']),
        f"{GENMAP_RESOURCE_DIR}/DG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/DG_genMap.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_genMap.csv"
    output:
        f"{FIGURES_DIR}/circos/refs_vs_genMap/data/DG_genMap_karyotype.txt",
        expand(f"{FIGURES_DIR}/circos/refs_vs_genMap/data/DG_{{match}}_genMap_links.txt", match = ['match', 'nomatch']),
        f"{FIGURES_DIR}/circos/refs_vs_genMap/data/DG_genMap_markerPositions.txt",
        f"{FIGURES_DIR}/circos/refs_vs_genMap/data/SG_genMap_karyotype.txt",
        expand(f"{FIGURES_DIR}/circos/refs_vs_genMap/data/SG_{{match}}_genMap_links.txt", match = ['match', 'nomatch']),
        f"{FIGURES_DIR}/circos/refs_vs_genMap/data/SG_genMap_markerPositions.txt"
    conda: '../envs/circos.yaml'
    script:
        "../scripts/r/generate_circos_genMap_links.R"

rule utm_TrRFive_genMap_circos:
    """
    Figure 1 and S2: Generate Circos plot with Linkage Markers aligned to both references, 
    with lines colored based on whether they align to the correct LG. Done for DG and SG populations
    """
    input:
        expand(rules.create_karyotype_file_ref.output, ref=['utm', 'TrRv5']),
        rules.generate_circos_genMap_links.output
    output:
        f"{FIGURES_DIR}/circos/refs_vs_genMap/utm_TrRFive_genMap{{map_pop}}_circos.png"
    log: f"{LOG_DIR}/circos/utm_TrRFive_genMap{{map_pop}}_circos.log"
    container: 'library://james-s-santangelo/circos/circos:v0.69-9'
    params:
        conf = lambda w: '../config/utm_TrRv5_vs_genMapDG_circos.conf' if w.map_pop == 'DG' else '../config/utm_TrRv5_vs_genMapSG_circos.conf',
        outdir = f"{FIGURES_DIR}/circos/refs_vs_genMap",
        outfile = 'utm_TrRFive_genMap{map_pop}_circos.png'
    shell:
        """
        circos -conf {params.conf} \
            -png \
            -nosvg \
            -outputfile {params.outfile} \
            -outputdir {params.outdir} 2> {log}
        """


##############################
#### UTM_TREP_v1.0 CIRCOS #### 
##############################

#### GC CONTENT ####

rule bedtools_nuc:
    """
    Use Bedtools nuc to estimate GC% in windows across the UTM reference
    """
    input:
        fasta = rules.split_fasta_toChroms_andOrganelles.output.chroms,
        bed = expand(rules.bedtools_makewindows.output.bed, ref = 'utm')
    output:
        f"{FIGURES_DIR}/circos/utm/data/utm_windowed_gc_content.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedtools/utm_nuc.log"
    shell:
        """
        bedtools nuc -fi {input.fasta} -bed {input.bed} |\
            tail -n +1 | cut -f1,2,3,6 > {output} 2> {log}
        """

#### GENE AND REPEAT DENSITY ####

rule gff_genesOnly:
    """
    Filter GFF to only include gene features
    """
    input:
        rules.gff_sort_functional.output 
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/utm/utm_genes.gff"
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
    """
    Convert Genes or Repeat (i.e. from Repeat Masker) GFF to BED files
    """
    input:
        lambda w: rules.repeat_masker.output.gff if w.feat == 'repeat' else rules.gff_genesOnly.output  
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/utm/utm_{{feat}}.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/utm_{{feat}}_toBed.log"
    shell:
        """
        gff2bed < {input} > {output} 2> {log}
        """

rule chromLengths_toBed:
    """
    Convert chromosmal length text files to BED format with start and end coordinates
    """
    input:
        bed = expand(rules.get_chrom_lengths.output, ref = 'utm')
    output:
         f"{PROGRAM_RESOURCE_DIR}/circos/utm_chrom_lengths.bed"
    conda: '../envs/circos.yaml'
    shell:
        """
        awk -vOFS="\t" '{{ print $1, "0", $2; }}' {input.bed} | sort-bed - > {output}
        """

rule bedops:
    """
    Use Bedops to generated bedmap-compatible windows file. Windows are 500 kb with 100 kb step
    """
    input:
        bed = rules.chromLengths_toBed.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/circos/utm_bedops_windows.bed"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/bedops/utm_bedops.log"
    shell:
        """
        bedops \
            --chop 500000 \
            --stagger 100000 \
            -x {input.bed} > {output} 2> {log}
        """

rule bedmap_featureCount:
    """
    Use Bedmap to estimate gene and repeat density in windows across the genome
    """
    input:
        win = rules.bedops.output,
        feat_bed = rules.gffToBed.output
    output:
        f"{FIGURES_DIR}/circos/utm/data/utm_{{feat}}.txt"
    conda: '../envs/circos.yaml'
    log: f"{LOG_DIR}/feature_count/utm_{{feat}}_count.log"
    params:
        val = lambda wildcards: 'count' if wildcards.feat == 'gene' else 'bases-uniq-f'
    shell:
        """
        bedmap \
            --echo \
            --{params.val} \
            --delim '\t' \
            {input.win} {input.feat_bed} > {output} 2> {log}
        """

rule utm_circos:
    """
    Figure 2: Generate Circos plot of UTM haploid mapping reference with GC%. gene density, and repeat density tracks. Used as template.
    Picture in center added afterwards in Inkscape.
    """
    input:
        expand(rules.create_karyotype_file_ref.output, ref = ['utm', 'TrRv5']),
        rules.bedtools_nuc.output,
        expand(rules.bedmap_featureCount.output, feat=['repeat', 'gene'])
    output:
        plot = f"{FIGURES_DIR}/circos/utm/utm_circos.png",
        rev = f"{FIGURES_DIR}/circos/utm/data/utm_karyotype_reversed.txt"
    log: f"{LOG_DIR}/circos/utm_TrRFive_genMap_circos.log"
    container: 'library://james-s-santangelo/circos/circos:v0.69-9'
    params:
        conf = '../config/utm_circos.conf',
        outdir = f"{FIGURES_DIR}/circos/utm",
        outfile = 'utm_circos.png'
    shell:
        """
        tac {input[0]} > {output.rev}
        circos -conf {params.conf} \
            -png \
            -nosvg \
            -param max_points_per_track=50000 \
            -outputfile {params.outfile} \
            -outputdir {params.outdir} 2> {log}
        """
      
rule circos_done:
    input:
        expand(rules.utm_TrRFive_genMap_circos.output, map_pop = ['SG', 'DG']),
        rules.utm_circos.output
    output:
        f"{FIGURES_DIR}/circos/circos.done"
    shell:
        """
        touch {output}
        """
        






