# Rules to annotate the genome

#####################################
#### REPEAT MODELING AND MASKING ####
#####################################

rule build_repeat_modeler_db:
    input:
        get_haplotype_fasta
    output:
        multiext(f"{ANNOTATION_DIR}/repeat_modeler/{{hap}}_rmdb", '.nhr', '.nin', '.nnd', '.nni', '.nog', '.nsq', '.translation')
    conda: '../envs/annotation.yaml'
    #container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/build_repeat_modeler_db/{hap}_rmdb.log'
    params:
        out = f"{ANNOTATION_DIR}/repeat_modeler/{{hap}}_rmdb"
    shell:
        """
        BuildDatabase -name {params.out} {input} 2> {log}
        """

rule repeat_modeler:
    input:
        rules.build_repeat_modeler_db.output
    output:
        fasta = f"{ANNOTATION_DIR}/repeat_modeler/{{hap}}_rmdb-families.fa",
        stk = f"{ANNOTATION_DIR}/repeat_modeler/{{hap}}_rmdb-families.stk"
    threads: 32
    conda: '../envs/annotation.yaml'
    #container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/repeat_modeler/{hap}_rm.log'
    params:
        db_base = lambda wildcards, input: os.path.splitext(input[0])[0]
    shell:
        """
        RepeatModeler -database {params.db_base} \
            -pa {threads} \
            -recoverDir ./{wildcards.hap} \
            -LTRStruct &> {log}
        """

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

rule merge_repeat_databases:
    input:
        rm_db = rules.configure_repbase.output.rm,
        tr_db = rules.repeat_modeler.output.fasta
    output:
        f"{PROGRAM_RESOURCE_DIR}/Libraries/{{hap}}_rm_merged_db.fasta"
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/merge_repeat_databases/{hap}_merge_repeat_databases.log'
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
        fasta = get_haplotype_fasta
    output:
        fasta = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/{{hap}}_softMasked.fasta",
        cat = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/{{hap}}_repeatMasker.cat.gz",
        out = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/{{hap}}_repeatMasker.out",
        gff = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/{{hap}}_repeatMasker.gff",
        stats = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/{{hap}}_repeatMasker.tbl"
    threads: 32
    #conda: '../envs/annotation.yaml'
    container: 'docker://dfam/tetools:1.6'
    log: LOG_DIR + '/repeat_masker/{hap}_repeat_masker.log'
    params:
        outdir = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/"
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

rule annotation_done:
    input:
        expand(rules.repeat_masker.output, hap='occ1')
    output:
        f"{ANNOTATION_DIR}/annotation.done"
    shell:
        """
        touch {output}
        """
