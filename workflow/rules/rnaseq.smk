# Temporary Snakefile with rules to map RNAseq reads and maybe do some BRAKER
# Eventually merge with annotation.smk

rule fasterq_dump:
    output:
        R1 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fq"),
        R2 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fq")
    log: LOG_DIR + '/fastq_dump/{{acc}}_fastq_dump.log'
    conda: '../envs/rnaseq.yaml'
    params:
        outdir = f"{ANNOTATION_DIR}/rnaseq_reads"
    shell:
        """
        fasterq-dump --split-3 \
            --skip-technical \
            --outdir {params.outdir} \
            {wildcards.acc} 2> {log}
        """

rule rnaseq_done:
    input:
        expand(rules.fasterq_dump.output, acc=RNASEQ_ACCESSIONS) 
    output:
        f"{ANNOTATION_DIR}/rna_seq.done"
    shell:
        """
        touch {output}
        """
