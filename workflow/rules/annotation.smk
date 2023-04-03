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
        R1 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fastq"),
        R2 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fastq")
    log: LOG_DIR + '/fastq_dump/{acc}_fastq_dump.log'
    conda: '../envs/annotation.yaml'
    params:
        outdir = f"{ANNOTATION_DIR}/rnaseq_reads"
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = lambda wildcards, attempt: str(attempt * 2) + ":00:00" 
    shell:
        """
        cd {params.outdir}
        fasterq-dump --split-3 \
            -e {threads} \
            --skip-technical \
            {wildcards.acc} 2> {log}
        """

rule gzip_fastq:
    input:
        rules.fasterq_dump.output
    output:
        R1 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fastq.gz"),
        R2 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fastq.gz")
    log: LOG_DIR + '/gzip_fastq/{acc}_gzip.log'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = lambda wildcards, attempt: str(attempt * 2) + ":00:00"
    shell:
        """
        gzip {input} 2> {log} 
        """

rule build_star:
    input:
        masked_genome = rules.repeat_masker.output.fasta 
    output:
        temp(directory(f"{ANNOTATION_DIR}/star/star_build"))
    log: LOG_DIR + '/star/star_build.log'
    conda:'../envs/annotation.yaml'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 15000,
        time = "1:00:00"
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.masked_genome} \
            --runThreadN {threads} &> {log}
        """

rule align_star:
    input:
        star_build = rules.build_star.output,
        R1 = rules.gzip_fastq.output.R1,
        R2 = rules.gzip_fastq.output.R2
    output:
        star_align = temp(f"{ANNOTATION_DIR}/star/star_align/{{acc}}/{{acc}}_Aligned.sortedByCoord.out.bam") 
    log: LOG_DIR + '/star/{acc}_star_align.log'
    conda:'../envs/annotation.yaml'
    params:
        out = f"{ANNOTATION_DIR}/star/star_align/{{acc}}/{{acc}}_"
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000,
        time = lambda wildcards, attempt: str(attempt * 12) + ":00:00"
    shell:
        """
        STAR --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --genomeDir {input.star_build} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --alignIntronMax 10000 \
            --outFileNamePrefix {params.out} &> {log} 
        """

###############################
#### STRUCTURAL ANNOTATION ####
###############################

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

rule download_uniprot_fabaceae_db:
    output:
        f'{PROGRAM_RESOURCE_DIR}/uniprot/fabaceae_proteins.fasta'
    log: LOG_DIR + '/uniprot_download/fabaceae_proteins_dl.log'
    params:
        url = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28taxonomy_name%3Afabaceae%29"
    shell:
        """
        curl --output {output} '{params.url}' 2> {log}
        """

rule combine_protein_dbs:
    input:
        ortho = rules.viridiplantae_orthodb.output,
        fab = rules.download_uniprot_fabaceae_db.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/allProteins.fasta"
    shell:
        """
        cat {input.fab} {input.ortho} > {output}
        """

rule braker_protein:
    input:
        proteins = rules.combine_protein_dbs.output,
        masked_genome = rules.repeat_masker.output.fasta,
    output:
        directory(f"{ANNOTATION_DIR}/braker/proteins")
    log: LOG_DIR + '/braker/proteins.log'
    params:
        genemark = GENEMARK,
        prothint = PROTHINT,
    threads: 30
    container: 'docker://teambraker/braker3' 
    shell:
        """
        braker.pl --genome {input.masked_genome} \
            --prot_seq {input.proteins} \
            --softmasking \
            --useexisting \
            --GENEMARK_PATH {params.genemark} \
            --PROTHINT_PATH {params.prothint} \
            --threads {threads} \
            --workingdir {output} \
            --species "Trifolium repens prot" 2> {log}
        """ 

rule count_uniprot_seqs:
    input:
        rules.download_uniprot_fabaceae_db.output
    output:
        f"{ANNOTATION_DIR}/uniprotSeqs_byTaxon.txt"
    conda: '../envs/notebooks.yaml'
    log: LOG_DIR + '/notebooks/count_uniprot_seqs_processed.ipynb'
    notebook:
        "../notebooks/count_uniprot_seqs.py.ipynb"

rule merge_rnaseq_bams:
    input:
        Star_Bam = expand(rules.align_star.output, acc=RNASEQ_ACCESSIONS),
    output:
        bam = f"{ANNOTATION_DIR}/star/allBams_merged.bam",
        bai = f"{ANNOTATION_DIR}/star/allBams_merged.bam.bai"
    conda: '../envs/annotation.yaml'
    threads: 30
    log: LOG_DIR + '/merge_rnaseq_bams/bams_merge.log'
    shell:
        """
        ( samtools merge -@ {threads} -r -o {output.bam} {input} &&\
            samtools index {output.bam} ) 2> {log}
        """

rule braker_rnaseq:
    input:
        masked_genome = rules.repeat_masker.output.fasta,
        bam = rules.merge_rnaseq_bams.output.bam
    output:
        directory(f"{ANNOTATION_DIR}/braker/rnaseq")
    log: LOG_DIR + '/braker/rnaseq.log'
    params:
        outputdir = f"{ANNOTATION_DIR}/braker/rnaseq",
        genemark=GENEMARK
    threads: 30
    container: 'docker://teambraker/braker3' 
    shell:
        """
        braker.pl --genome {input.masked_genome} \
            --bam {input.bam} \
            --softmasking \
            --useexisting \
            --threads {threads} \
            --GENEMARK_PATH={params.genemark} \
            --workingdir {output} \
            --species "Trifolium repens rna" 2> {log} 
        """

rule tsebra_combine:
    input:
        rules.braker_protein.output,
        rules.braker_rnaseq.output,
        config = '../config/tsebra.cfg',
        rna_aug = f"{rules.braker_rnaseq.output[0]}/Augustus/augustus.hints.gtf",
        prot_aug = f"{rules.braker_protein.output[0]}/Augustus/augustus.hints.gtf",
        hints_rna= f"{rules.braker_rnaseq.output[0]}/hintsfile.gff",
        hints_prot = f"{rules.braker_protein.output[0]}/hintsfile.gff",
    output:
        braker_combined = f"{ANNOTATION_DIR}/braker/tsebra/braker_combined.gtf"
    log: LOG_DIR + '/braker/tsebra.log'
    container: 'docker://teambraker/braker3' 
    shell:
        """
        tsebra.py -g {input.rna_aug},{input.prot_aug} \
            -c {input.config} \
            -e {input.hints_rna},{input.hints_prot} \
            -o {output} 2> {log} 
        """ 

rule rename_tsebra_gtf:
    input:
        rules.tsebra_combine.output
    output:
        gtf = f"{ANNOTATION_DIR}/braker/tsebra/braker_combined_renamed.gtf",
        tab = f"{ANNOTATION_DIR}/braker/tsebra/tsebra_rename_translationTab.txt"
    log: LOG_DIR + '/braker/rename_gtf.log'
    container: 'docker://teambraker/braker3'
    shell:
        """
        rename_gtf.py --gtf {input} \
            --prefix TrR_v6 \
            --translation_tab {output.tab} \
            --out {output.gtf} 2> {log} 
        """

rule addUTRs:
    input:
        masked_genome = rules.repeat_masker.output.fasta,
        rna_bam = rules.merge_rnaseq_bams.output.bam,
        hints = rules.rename_tsebra_gtf.output.gtf
    output:
        f"{ANNOTATION_DIR}/braker/gushr/TrR_v6_brakerCombined_tsebraRenamed_withUTRs.gtf"
    log: LOG_DIR + '/braker/gushr_addUTRs.log'
    params:
        tmp = f"{config['results_prefix']}/tmp/",
        out = f"{ANNOTATION_DIR}/braker/gushr/TrR_v6_brakerCombined_tsebraRenamed_withUTRs",
        gemoma = f"{config['results_prefix']}/tmp/GeMoMa-1.6.2.jar"
    container: 'library://james-s-santangelo/gushr/gushr:v1.0.0'
    threads: 10
    shell:
        """
        mkdir -p {params.tmp}
        cp /opt/GUSHR/GeMoMa-1.6.2.jar {params.gemoma}
        gushr.py  -t {input.hints} \
            -b {input.rna_bam} \
            -g {input.masked_genome} \
            -o {params.out} \
            -q 1 \
            -c {threads} \
            -d {params.tmp} \
            --GeMoMaJar {params.gemoma}
        """

rule clean_gushr_gtf:
    input:
        rules.addUTRs.output
    output:
        f"{ANNOTATION_DIR}/braker/gushr/TrR_v6_brakerCombined_tsebraRenamed_withUTRs_cleaned.gtf"
    run:
        with open(input[0], 'r') as fin:
            with open(output[0], 'w') as fout:
                lines = fin.readlines()
                for l in lines:
                    sl = l.split('\t')
                    ssl = list(map(str.strip, sl))
                    ssl_3pr = ["three_prime_utr" if x=="3'-UTR" else x for x in ssl]
                    ssl_5pr = ["five_prime_utr" if x=="5'-UTR" else x for x in ssl_3pr]
                    nssl = "\t".join(ssl_5pr)
                    fout.write(f"{nssl}\n") 

###############################
#### CLEAN & TRANSFORM GTF ####
###############################

rule gtf_to_gff:
    input:
        rules.clean_gushr_gtf.output
    output:
        f"{ANNOTATION_DIR}/TrR_v6_structural.gff"
    container: 'docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0'
    log: LOG_DIR + '/gtf_to_gff/gtf_to_gff.log'
    shell:
        """
        agat_convert_sp_gxf2gxf.pl --gtf {input} \
            --output {output} &> {log}
        """

rule gff_sort:
    input:
        rules.gtf_to_gff.output
    output:
        f"{ANNOTATION_DIR}/TrR_v6_structural_sorted.gff"
    log: LOG_DIR + '/gff_sort/gff_sort.log'
    conda: '../envs/annotation.yaml'
    shell:
        """
        gff3sort.pl {input} > {output} 2> {log}
        """

rule get_proteins:
    input:
        gff = rules.gff_sort.output,
        ref = rules.repeat_masker.output.fasta
    output:
        f"{ANNOTATION_DIR}/TrR_v6_proteins.fasta"
    log: LOG_DIR + '/get_proteins/get_proteins.log'
    container: 'docker://teambraker/braker3'
    shell:
        """
        gffread -E -y {output} -g {input.ref} {input.gff} 2> {log}
        """

###############################
#### FUNCTIONAL ANNOTATION ####
###############################

rule run_interproscan:
    input:
        data = IPRSCAN_DATA,
        prot = rules.get_proteins.output
    output:
        gff =  f"{ANNOTATION_DIR}/interproscan/TrR_v6_interproscan.gff3",
        xml =  f"{ANNOTATION_DIR}/interproscan/TrR_v6_interproscan.xml"
    log: LOG_DIR + '/interproscan/run_interproscan.log'
    threads: 32
    params:
        out_base =  f"{ANNOTATION_DIR}/interproscan/TrR_v6_interproscan"
    container: 'library://james-s-santangelo/interproscan/interproscan:5.61-93.0' 
    shell:
        """
        interproscan.sh -i {input.prot} \
            -b {params.out_base} \
            -f xml,gff3 \
            -goterms \
            --pathways \
            --seqtype p \
            --cpu {threads} \
            --verbose &> {log} 
        """ 

rule funannotate_setup:
    output:
        directory(f"{ANNOTATION_DIR}/funannotate/fun_db")
    log: LOG_DIR + '/funannotate/funannotate_setup.log'
    container: '/home/santang3/singularity_containers/funannotate.sif'
    shell:
        """
        funannotate setup --database {output} -b embryophyta --force 2> {log}
        """

rule feature_recode:
    input:
        rules.gff_sort.output
    output:
        f"{ANNOTATION_DIR}/TrR_v6_structural_sorted_featureRecode.gff"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\t"}} {{gsub(/transcript/, "mRNA", $3)}} 1' {input} > {output}
        """


rule dl_eggnog_db:
    output:
        directory(f"{ANNOTATION_DIR}/eggnog/eggnog_db")
    conda: '../envs/eggnog.yaml'
    shell:
        """
        mkdir -p {output}
        download_eggnog_data.py -y --data_dir {output} 
        """

rule run_eggnog_mapper:
    input:
        prot = rules.get_proteins.output,
        db = rules.dl_eggnog_db.output
    output:
        annot = f"{ANNOTATION_DIR}/eggnog/TrR_v6.emapper.annotations"
    log: LOG_DIR + '/eggnog/aggnog_mapper.log'
    threads: 24
    conda: '../envs/eggnog.yaml'
    params:
        out = f"{ANNOTATION_DIR}/eggnog/TrR_v6"
    shell:
        """
        emapper.py -i {input.prot} \
            --data_dir {input.db} \
            --cpu {threads} \
            --itype proteins \
            -o {params.out} \
            --override &> {log}

        """

rule funannotate_annotate:
    input:
        enm = rules.run_eggnog_mapper.output.annot,
        gff = rules.feature_recode.output,
        ref = rules.repeat_masker.output.fasta,
        iprs = rules.run_interproscan.output.xml,
        db = rules.funannotate_setup.output 
    output:
        directory(f"{ANNOTATION_DIR}/funannotate/annotations")
    log: LOG_DIR + '/funannotate/funannotate_annotate.log'
    container: '/home/santang3/singularity_containers/funannotate.sif'
    threads: 32
    params:
        sbt = NCBI_TEMPLATE,
        locus_tag = 'P8452'
    shell:
        """
        funannotate annotate \
            --sbt {params.sbt} \
            --gff {input.gff} \
            --fasta {input.ref} \
            --species "Trifolium repens" \
            --out {output} \
            --iprscan {input.iprs} \
            --eggnog {input.enm} \
            --force \
            --rename {params.locus_tag} \
            --database {input.db} \
            --busco_db embryophyta \
            --cpus {threads} 2> {log}
        """

rule annotation_done:
    input:
        expand(rules.funannotate_annotate.output)
    output:
        f"{ANNOTATION_DIR}/annotation.done"
    shell:
        """
        touch {output}
        """
