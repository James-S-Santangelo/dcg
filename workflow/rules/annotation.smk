# Rules to annotate the genome

#####################################
#### REPEAT MODELING AND MASKING ####
#####################################

rule configure_repbase:
    input:
        REPBASE
    output:
        libdir = directory(f"{PROGRAM_RESOURCE_DIR}/Libraries"),
        dfam = f"{PROGRAM_RESOURCE_DIR}/Libraries/Dfam.h5",
        rm = f"{PROGRAM_RESOURCE_DIR}/Libraries/RepeatMaskerLib.h5"
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/configure_repbase/configure_repbase.log'
    params:
        untar_dir = f"{PROGRAM_RESOURCE_DIR}"
    shell:
        """
        ( tar -xzf {input} -C {params.untar_dir} &&
        cp -r /opt/RepeatMasker/Libraries/* {output.libdir} &&
        ln -sf {output.dfam} {output.rm} &&
        addRepBase.pl --libdir {output.libdir} ) &> {log}
        """

rule build_repeat_modeler_db:
    input:
        rules.create_reference_assemblies.output
    output:
        multiext(f"{ANNOTATION_DIR}/repeat_modeler/rmdb", '.nhr', '.nin', '.nnd', '.nni', '.nog', '.nsq', '.translation')
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/build_repeat_modeler_db/rmdb.log'
    params:
        out = f"{ANNOTATION_DIR}/repeat_modeler/rmdb"
    shell:
        """
        BuildDatabase -name {params.out} {input} &> {log}
        """

rule repeat_modeler:
    input:
        rules.build_repeat_modeler_db.output,
        rules.configure_repbase.output
    output:
        fasta = f"{ANNOTATION_DIR}/repeat_modeler/rmdb-families.fa",
        stk = f"{ANNOTATION_DIR}/repeat_modeler/rmdb-families.stk"
    threads: 48
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/repeat_modeler/rm.log'
    params:
        db_base = lambda wildcards, input: os.path.splitext(input[0])[0]
    shell:
        """
        RepeatModeler -database {params.db_base} \
            -pa {threads} \
            -LTRStruct &> {log}
        """

rule merge_repeat_databases:
    input:
        rm_db = rules.configure_repbase.output.rm,
        tr_db = rules.repeat_modeler.output.fasta
    output:
        f"{PROGRAM_RESOURCE_DIR}/Libraries/rm_merged_db.fasta"
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/merge_repeat_databases/merge_repeat_databases.log'
    params:
        rm_db_fasta = f"{PROGRAM_RESOURCE_DIR}/Libraries/rm_db.fasta"
    shell:
        """
        ( famdb.py -i {input.rm_db} families \
            --format fasta_name \
            --ancestors --descendants --include-class-in-name \
            'viridiplantae' > {params.rm_db_fasta} &&

        cat {params.rm_db_fasta} {input.tr_db} > {output} ) 2> {log}
        """

rule repeat_masker:
    input:
        lib = rules.merge_repeat_databases.output,
        fasta = rules.create_reference_assemblies.output
    output:
        fasta = f"{ANNOTATION_DIR}/repeat_masker/TrR_v6_haploid_reference_softMasked.fasta",
        cat = f"{ANNOTATION_DIR}/repeat_masker/TrR_v6_haploid_reference_repeatMasker.cat.gz",
        out = f"{ANNOTATION_DIR}/repeat_masker/TrR_v6_haploid_reference_repeatMasker.out",
        gff = f"{ANNOTATION_DIR}/repeat_masker/TrR_v6_haploid_reference_repeatMasker.gff",
        stats = f"{ANNOTATION_DIR}/repeat_masker/TrR_v6_haploid_reference_repeatMasker.tbl"
    threads: 48
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/repeat_masker/repeat_masker.log'
    params:
        outdir = f"{ANNOTATION_DIR}/repeat_masker/"
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

###########################################
#### DOWNLOAD AND MAP ALL RNASEQ READS ####
###########################################

rule prefetch:
    output:
        temp(directory(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}"))
    log: LOG_DIR + '/prefetch/{acc}.log'
    conda: '../envs/annotation.yaml'
    params:
        outdir = f"{ANNOTATION_DIR}/rnaseq_reads"
    shell:
        """
        prefetch {wildcards.acc} -O {params.outdir} 2> {log}
        """

rule fasterq_dump:
    input:
        expand(rules.prefetch.output, acc=RNASEQ_ACCESSIONS)
    output:
        R1 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fq"),
        R2 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fq")
    log: LOG_DIR + '/fastq_dump/{acc}_fastq_dump.log'
    conda: '../envs/annotation.yaml'
    params:
        outdir = f"{ANNOTATION_DIR}/rnaseq_reads"
    threads: 6
    resources:
        mem_mb = lambda attempt, wildcards: attempt * 2000,
        time = "1:00:00"
    shell:
        """
        cd {params.outdir}
        fasterq-dump --split-3 \
            -e {threads} \
            --skip-technical \
            {wildcards.acc} 2> {log}
        """

rule build_star:
    input:
        masked_genome = rules.repeat_masker.output.fasta 
    output:
        temp(directory(f"{ANNOTATION_DIR}/star/star_build"))
    log: LOG_DIR + '/star/star_build.log'
    conda:'../envs/annotation.yaml'
    threads: 16
    resources:
        mem_mb = lambda attempt, wildcards: attempt * 2000,
        time = "1:00:00"
    shell:
        """
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFile {input.masked_genome} \
            --runThreadN {threads} 2> {log}
        """

rule align_star:
    input:
        unpack(get_star_align_input_files)
    output:
        star_align = f"{ANNOTATION_DIR}/star/star_align/{{acc}}_sorted.bam" 
    log: LOG_DIR + '/star/{acc}_star_align.log'
    conda:'../envs/annotation.yaml'
    params:
        out = f"{ANNOTATION_DIR}/star/star_align/{{acc}}_sorted"
    threads: 6
    resources:
        mem_mb = lambda attempt, wildcards: attempt * 10000,
        time = "6:00:00"
    shell:
        """
        STAR --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --genomeDir {input.star_build} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.out} 2> {log} 
        """

################
#### BRAKER ####
################

rule viridiplantae_orthodb:
    output:
        Plant_ProrteinDB = f"{PROGRAM_RESOURCE_DIR}/orthodb/Viridiplantae_protein.fasta"
    log: LOG_DIR + "/orthodb/ortho.log"
    params:
        outdir = f"{PROGRAM_RESOURCE_DIR}/orthodb"
    shell:
        """
        ( wget --no-check-certificate https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz && 
        tar -zxf odb10_plants_fasta.tar.gz -C {params.outdir} &&
        cat {params.outdir}/plants/Rawdata/* > {output} ) 2> {log}
        """

rule braker_protein:
    input:
        proteins = rules.viridiplantae_orthodb.output,
        masked_genome = rules.repeat_masker.output.fasta,
    output:
        hints_protein = f"{ANNOTATION_DIR}/braker/proteins/hintsfile.gff",
        aug_hint_protein = f"{ANNOTATION_DIR}/braker/proteins/augustus.hints.gtf"
    log: LOG_DIR + '/braker/braker_proteins.log'
    params:
        outputdir = f"{ANNOTATION_DIR}/braker/proteins",
        aug_config = "../resources/augustus_config"
    threads: 20
    container: 'library://james-s-santangelo/braker/braker:2.1.6'
    shell:
        """
        export AUGUSTUS_CONFIG_PATH={params.aug_config}
        braker.pl --genome {input.masked_genome} \
            --prot_seq {input.proteins} \
            --epmode \
            --softmasking \
            --cores {threads} \
            --workingdir {params.outputdir} \
            --species "Trifolium repens" 2> {log}
        """ 

rule braker_rnaseq:
    input:
        masked_genome = rules.repeat_masker.output.fasta,
        Star_Bam = expand(rules.align_star.output, acc=ALL_RNASEQ_SAMPLES, hap='occ1')
    output:
        hints_rna = f"{ANNOTATION_DIR}/braker/braker_rnaseq/hintsfile.gff", 
        aug_hint_rna = f"{ANNOTATION_DIR}/braker/braker_rnaseq/augustus.hints.gtf"
    log: LOG_DIR + '/braker/braker_annotate.log'
    params:
        outputdir = f"{ANNOTATION_DIR}/braker/braker_rnaseq",
        aug_config = "../resources/augustus_config"
    threads: 20
    container: 'library://james-s-santangelo/braker/braker:2.1.6'
    shell:
        """
        export AUGUSTUS_CONFIG_PATH={params.aug_config}
        braker.pl --genome {input.masked_genome} \
            --bam {input.Star_Bam} \
            --softmasking \
            --cores {threads} \
            --workingdir {params.outputdir} \
            --species "Trifolium repens" 2> {log} 
        """

rule tsebra_combine:
    input:
        rna_aug = rules.braker_rnaseq.output.aug_hint_rna,
        protein_aug = rules.braker_protein.output.aug_hint_protein,
        hints_rna = rules.braker_rnaseq.output.hints_rna,
        hints_protein = rules.braker_protein.output.hints_protein
    output:
        braker_combined = f"{ANNOTATION_DIR}/braker/tsebra/braker_combined.gtf"
    log: LOG_DIR + '/braker/tsebra.log'
    container: 'library://james-s-santangelo/braker/braker:2.1.6'
    shell:
        """
        tesbra.py -g {input.rna_aug},{input.protein_aug} \
            -c default.cfg \
            -e {input.hints_rna},{input.hints_protein} \
            -o {output} 2> {log} 
        """ 

rule annotation_done:
    input:
        expand(rules.fasterq_dump.output, acc=RNASEQ_ACCESSIONS),
        rules.tsebra_combine.output
    output:
        f"{ANNOTATION_DIR}/annotation.done"
    shell:
        """
        touch {output}
        """
