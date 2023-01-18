# Rules for finding and working with HCN loci in the new Haplotype-resolved assembly

rule blast_hcn_loci:
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = HCN_LOCI_FASTA
    output:
        f"{ORGANELLE_HCN_DIR}/hcn_loci/{{hap}}_hcnLoci_blast.txt"
    conda: '../envs/organelles_hcn.yaml'
    log: LOG_DIR + '/organelle_hcn/{hap}_hcnLoci_blast.log'
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

rule blast_organelle_seqs:
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = ORGANELLE_SEQS
    output:
        f"{ORGANELLE_HCN_DIR}/organelles/{{hap}}_organelles_blast.txt"
    conda: '../envs/organelles_hcn.yaml'
    log: LOG_DIR + '/organelle_hcn/{hap}_organelles_blast.log'
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

rule organelles_hcn_done:
    input:
        expand(rules.blast_hcn_loci.output, hap=HAPS),
        expand(rules.blast_organelle_seqs.output, hap=HAPS)
    output:
        f"{ORGANELLE_HCN_DIR}/organelles_hcn.done"
    shell:
        """
        touch {output}
        """
