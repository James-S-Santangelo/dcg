# Rules for finding and working with HCN loci in the new Haplotype-resolved assembly

rule blast_hcn_loci:
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = HCN_LOCI_FASTA
    output:
        f"{BLAST_DIR}/hcn_loci/{{hap}}_hcnLoci_blast.txt"
    conda: '../envs/blast.yaml'
    log: LOG_DIR + '/blast/{hap}_hcnLoci_blast.log'
    threads: 4
    params:
        outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp'",
        db_base = lambda wildcards, input: os.path.splitext(input.db[0])[0]
    shell:
        """
        blastn -db {params.db_base} \
            -query {input.query} \
            -out {output} \
            -outfmt {params.outfmt} \
            -num_threads {threads} \
            -evalue 1e-10 \
            -max_hsps 5 \
            -max_target_seqs 5 2> {log}  
        """

rule blast_rbcl_coi:
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = RBCL_COI
    output:
        f"{BLAST_DIR}/rbcl_coi/{{hap}}_rbcl_coi_blast.txt"
    conda: '../envs/blast.yaml'
    log: LOG_DIR + '/blast/{hap}_rbcl_coi_blast.log'
    threads: 4
    params:
        outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp'",
        db_base = lambda wildcards, input: os.path.splitext(input.db[0])[0]
    shell:
        """
        blastn -db {params.db_base} \
            -query {input.query} \
            -out {output} \
            -outfmt {params.outfmt} \
            -num_threads {threads} \
            -evalue 1e-10 \
            -max_hsps 5 \
            -max_target_seqs 5 2> {log}  
        """

rule blast_done:
    input:
        expand(rules.blast_hcn_loci.output, hap=HAPS),
        expand(rules.blast_rbcl_coi.output, hap=HAPS)
    output:
        f"{BLAST_DIR}/blast.done"
    shell:
        """
        touch {output}
        """
