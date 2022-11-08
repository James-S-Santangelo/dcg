# Rules for fitting and interpolating recombination maps for new Dovetail assemblies

rule blast_markers:
    input:
        markers = ancient(f'{GENMAP_RESOURCE_DIR}/{{map_pop}}_tags.fa'),
        db_flag = rules.makeblastdb_fromHaplotypeFasta.output,
        db = rules.makeblastdb_fromHaplotypeFasta.output 
    output:
        f'{GENMAP_RESULTS_DIR}/{{hap}}_{{map_pop}}_marker_blast.txt'
    conda: '../envs/recombination_map.yaml'
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

rule gen_map_done:
    input:
        expand(rules.blast_markers.output, hap=HAPS, map_pop=['SG','DG'])
    output:
        f"{GENMAP_RESULTS_DIR}/gen_map.done"
    shell:
        """
        touch {output}
        """
