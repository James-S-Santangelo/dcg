# Rules for finding and working with HCN loci in the new Haplotype-resolved assembly

rule makeblastdb_fromHaplotypeFasta:
    """
    Creates BLAST Database from revised haplotype sequences.
    """
    input:
        f"{REVISED_HAP_DIR}/{{hap}}_revised.fasta"
    output:
        multiext(f'{BLAST_DIR}/blastDBs/{{hap}}/{{hap}}.fasta', '.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto') 
    conda: '../envs/blast.yaml'
    log: LOG_DIR + '/makeblastdb/makeblastdb_{hap}.log'
    params:
        outfile = f'{BLAST_DIR}/blastDBs/{{hap}}/{{hap}}.fasta'
    shell:
        """
        makeblastdb -in {input} \
            -dbtype nucl \
            -parse_seqids \
            -logfile {log} \
            -out {params.outfile}
        """

rule blast_markers:
    """
    BLASTs linkage markers from Olsen et al. (2022) F2 mapping populations against the revised haplptypes. Used to determine linkage groups and chromosomes. 
    """
    input:
        markers = ancient(f'{GENMAP_RESOURCE_DIR}/{{map_pop}}_tags.fa'),
        db_flag = rules.makeblastdb_fromHaplotypeFasta.output,
        db = rules.makeblastdb_fromHaplotypeFasta.output 
    output:
        f'{BLAST_DIR}/genMap/{{hap}}_{{map_pop}}_marker_blast.txt'
    conda: '../envs/blast.yaml'
    log: LOG_DIR + '/blast_markers/{hap}_{map_pop}_marker_blast.log'
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

rule blast_hcn_loci:
    """
    BLAST known sequences for HCN loci against the revised haplptypes to identify which haplotypes contain functional copies of Ac and Li.
    """
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = HCN_LOCI_FASTA
    output:
        f"{BLAST_DIR}/hcn_loci/{{hap}}_hcnLoci_blast.txt"
    conda: '../envs/blast.yaml'
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
    """
    BLAST known mitochondrial and chloroplast sequences against the revised haplotypes to identify which scaffolds correspond to the organelles in each haplotype FASTA.
    """
    input:
        db = rules.makeblastdb_fromHaplotypeFasta.output,
        query = ORGANELLE_SEQS
    output:
        f"{BLAST_DIR}/organelles/{{hap}}_organelles_blast.txt"
    conda: '../envs/blast.yaml'
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

rule blast_done:
    input:
        expand(rules.blast_hcn_loci.output, hap=HAPS),
        expand(rules.blast_organelle_seqs.output, hap=HAPS),
        expand(rules.blast_markers.output, hap=HAPS, map_pop=['SG','DG'])
    output:
        f"{BLAST_DIR}/blast.done"
    shell:
        """
        touch {output}
        """
