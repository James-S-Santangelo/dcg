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
        rules.create_reference_assemblies.output[0]
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
        fasta = rules.create_reference_assemblies.output[0]
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
        prot_aug = f"{ANNOTATION_DIR}/braker/proteins/Augustus/augustus.hints.gtf",
        hints_prot = f"{ANNOTATION_DIR}/braker/proteins/hintsfile.gff"
    log: LOG_DIR + '/braker/proteins.log'
    params:
        outputdir = f"{ANNOTATION_DIR}/braker/rnaseq",
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
            --workingdir {params.outputdir} \
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
        rna_aug = f"{ANNOTATION_DIR}/braker/rnaseq/Augustus/augustus.hints.gtf",
        hints_rna = f"{ANNOTATION_DIR}/braker/rnaseq/hintsfile.gff"
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
            --workingdir {params.outputdir} \
            --species "Trifolium repens rna" 2> {log} 
        """

rule tsebra_combine:
    input:
        rules.braker_protein.output,
        rules.braker_rnaseq.output,
        rna_aug = rules.braker_rnaseq.output.rna_aug, 
        prot_aug = rules.braker_protein.output.prot_aug, 
        hints_rna= rules.braker_rnaseq.output.hints_rna,
        hints_prot = rules.braker_protein.output.hints_prot 
    output:
        braker_combined = f"{ANNOTATION_DIR}/braker/tsebra/braker_combined.gtf"
    log: LOG_DIR + '/braker/tsebra.log'
    container: 'docker://teambraker/braker3'
    params:
        config = '../config/tsebra.cfg',
    shell:
        """
        tsebra.py -g {input.rna_aug},{input.prot_aug} \
            -c {params.config} \
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
            --prefix ACLI19 \
            --translation_tab {output.tab} \
            --out {output.gtf} 2> {log} 
        """

###############################
#### CLEAN & TRANSFORM GTF ####
###############################

rule remove_features_and_organelles:
    input:
        rules.rename_tsebra_gtf.output
    output:
        f"{ANNOTATION_DIR}/cleaned/braker_combined_CDSonly_noOrgs.gtf"
    run:
        organelles = ['Mitochondria', 'Plastid']
        features = ['exon', 'intron', 'gene', 'transcript']
        with open(input[0], 'r') as fin:
            with open(output[0], 'w') as fout:
                lines = fin.readlines()
                for line in lines:
                    sline = line.split('\t')
                    if sline[0] in organelles or sline[2] in features:
                        pass
                    else:
                        fout.write(line)

rule gtf_to_gff:
    input:
        rules.remove_features_and_organelles.output
    output:
        f"{ANNOTATION_DIR}/cleand/UTM_Trep_v1.0_structural.gff"
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
        f"{ANNOTATION_DIR}/cleaned/UTM_Trep_v1.0_structural_sorted.gff"
    log: LOG_DIR + '/gff_sort/gff_sort.log'
    conda: '../envs/annotation.yaml'
    shell:
        """
        gt gff3 -sort -tidy -retainids {input} > {output} 2> {log}
        """

rule reformat_gff:
    input:
        rules.gff_sort.output
    output:
        f"{ANNOTATION_DIR}/cleaned/UTM_Trep_v1.0_structural_sorted_reformated.gff"
    run:
        with open(input[0], 'r') as fin:
            with open(output[0], 'w') as fout:
                lines = fin.readlines()
                for line in lines:
                    if not line.startswith('#'):
                        sline = line.split('\t')
                        feature = sline[2]
                        if feature == 'gene':
                            # Remove transcript_id attribute from gene features
                            sline[8] = re.sub(r'(;transcript_id.*$)', '', sline[8])
                        elif feature == 'mRNA':
                            # Make sure transcript_id annotations for isoforms are correct (i.e., .t2, .t3, etc.)
                            # Use ID attribute since transcript_id attributes are incorrectly incremented and will be replaced
                            id_pattern = r"(?<=ID=)(.*)(?=;Parent)"
                            ID = re.search(id_pattern, sline[8]).group(1)
                            if ID.endswith('.t1'):
                                # First isoforms are fine
                                pass
                            else:
                                # Alternative isoforms need transcript_id replaced with ID
                                sline[8] = re.sub(r'(?<=;transcript_id=)(.*$)', ID, sline[8])
                        fout.write('\t'.join(sline))
                    else:
                        fout.write(line) 


rule get_proteins:
    input:
        gff = rules.reformat_gff.output,
        ref = rules.repeat_masker.output.fasta
    output:
        f"{ANNOTATION_DIR}/UTM_Trep_v1.0_proteins.fasta"
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
    threads: 24
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

rule split_fasta_toChroms_andOrganelles:
    input:
       rules.repeat_masker.output.fasta
    output:
        chroms = f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_chromsOnly.fasta", 
        mito = f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_mitochondria.fasta", 
        cp = f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_chloroplast.fasta"
    conda: '../envs/annotation.yaml'
    params:
        chroms = 'Chr01_Occ Chr01_Pall Chr02_Occ Chr02_Pall Chr03_Occ Chr03_Pall Chr04_Occ Chr04_Pall Chr05_Occ Chr05_Pall Chr06_Occ Chr06_Pall Chr07_Occ Chr07_Pall Chr08_Occ Chr08_Pall'
    shell:
        """
        samtools faidx {input} {params.chroms} > {output.chroms}
        samtools faidx {input} Mitochondria > {output.mito}
        samtools faidx {input} Plastid > {output.cp}
        """

rule funannotate_annotate:
    input:
        enm = rules.run_eggnog_mapper.output.annot,
        gff = rules.reformat_gff.output, 
        ref = rules.split_fasta_toChroms_andOrganelles.output.chroms,
        iprs = rules.run_interproscan.output.xml,
        db = rules.funannotate_setup.output 
    output:
        directory(f"{ANNOTATION_DIR}/funannotate/annotations"),
        fasta = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Trifolium_repens.scaffolds.fa", 
        prot = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Trifolium_repens.proteins.fa", 
        gff3 = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Trifolium_repens.gff3", 
        agp = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Trifolium_repens.agp", 
        gbk = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Trifolium_repens.gbk", 
        tbl = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Trifolium_repens.tbl"
    log: LOG_DIR + '/funannotate/funannotate_annotate.log'
    container: 'docker://nextgenusfs/funannotate:v1.8.15'
    threads: 48 
    params:
        sbt = NCBI_TEMPLATE,
        outdir = f"{ANNOTATION_DIR}/funannotate/annotations"
    shell:
        """
        mkdir {params.outdir}
        funannotate annotate \
            --sbt {params.sbt} \
            --gff {input.gff} \
            --fasta {input.ref} \
            --species "Trifolium repens" \
            --out {params.outdir} \
            --iprscan {input.iprs} \
            --eggnog {input.enm} \
            --force \
            --database {input.db} \
            --busco_db embryophyta \
            --cpus {threads} 2> {log}
        """

rule gff_sort_functional:
    input:
        rules.funannotate_annotate.output.gff3
    output:
        f"{ANNOTATION_DIR}/UTM_Trep_v1.0_functional_sorted.gff"
    log: LOG_DIR + '/gff_sort/gff_sort_functional.log'
    conda: '../envs/annotation.yaml'
    shell:
        """
        gt gff3 -sort -tidy -retainids {input} > {output} 2> {log}
        """

rule add_locus_tag:
    input:
        rules.gff_sort_functional.output
    output:
        f"{ANNOTATION_DIR}/UTM_Trep_v1.0_functional_sorted_wLocusTag.gff"
    params:
        locus_tag = 'P8452'
    run:
        with open(input[0], 'r') as fin:
            with open(output[0], 'w') as fout:
                lines = fin.readlines()
                total_genes = sum([l.split('\t')[2] == 'gene' for l in lines if not l.startswith('#')])
                gene_num = 1
                for line in lines:
                    if not line.startswith('#'):
                        sline = line.split('\t')
                        feature = sline[2]
                        if feature == 'gene':
                            # Remove transcript_id attribute from gene features
                            sline[8] = re.sub(r'(;transcript_id.*$)', '', sline[8])
                                                
                            # Assign locus tag to attributes of gene features
                            locus_tag = f'{params.locus_tag}_' + str(int(gene_num)).zfill(len(str(total_genes)))
                            sline[8] = f"{sline[8].strip()};locus_tag={locus_tag}\n"
                            gene_num += 1
                        else:
                            pass
                        fout.write('\t'.join(sline))
                    else:
                        fout.write(line) 

rule tableToAsn_haploid:
    input:
        gff = rules.add_locus_tag.output,
        sbt = NCBI_TEMPLATE,
        ref = rules.split_fasta_toChroms_andOrganelles.output.chroms
    output:
        directory(f"{NCBI_DIR}/table2asn_haploid")
    log: f"{LOG_DIR}/tableToAsn/tableToAsn_haploid.log"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz
        gunzip linux64.table2asn.gz
        mv linux64.table2asn table2asn
        chmod +x table2asn

        ./table2asn -i {input.ref} \
            -f {input.gff} \
            -t {input.sbt} \
            -o {output} \
            -Z -V vb -euk -c f \
            -gaps-min 10 \
            -l proximity-ligation \
            -gaps-unknown 100 \
            -j "[organism=Trifolium repens]" 2> {log}
        """

rule annotation_done:
    input:
        expand(rules.tableToAsn_haploid.output)
    output:
        f"{ANNOTATION_DIR}/annotation.done"
    shell:
        """
        touch {output}
        """
