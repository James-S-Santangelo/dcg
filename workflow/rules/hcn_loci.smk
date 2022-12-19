# Rules for finding and working with HCN loci in the new Haplotype-resolved assembly

rule blast_hcn_loci:
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = HCN_LOCI_FASTA
    output:
        f"{HCN_LOCI_DIR}/blast/{{hap}}_blastResults.txt"
    conda: '../envs/hcn_loci.yaml'
    log: LOG_DIR + '/blast_hcn_loci/{hap}_blast.log'
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

rule hcn_loci_done:
    input:
        expand(rules.blast_hcn_loci.output, hap=HAPS)
    output:
        f"{HCN_LOCI_DIR}/hcn_loci.done"
    shell:
        """
        touch {output}
        """
