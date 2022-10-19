# Rules to work with the raw haplotype fastas returned by Dovetail

rule makeblastdb_fromHaplotypeFasta:
    input:
        get_haplotype_fasta
    output:
        multiext(f'{REF_DIR}/blastDBs/{{hap}}/{{hap}}.fasta', '.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto') 
    conda: '../envs/ref.yaml'
    log: LOG_DIR + '/makeblastdb/makeblastdb_{hap}.log'
    params:
        outfile = f'{REF_DIR}/blastDBs/{{hap}}/{{hap}}.fasta'
    shell:
        """
        makeblastdb -in {input} \
            -dbtype nucl \
            -parse_seqids \
            -logfile {log} \
            -out {params.outfile}
        """

rule ref_done:
    input:
        expand(rules.makeblastdb_fromHaplotypeFasta.output, hap=HAPS)
    output:
        f"{REF_DIR}/ref.done"
    shell:
        """
        touch {output}
        """

