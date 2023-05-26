rule keep_longest_isoform:
    """
    Uses GFFutil module to iterate through features and keep only the 
    longest isoform for each gene
    """
    input:
        gff = rules.gff_sort_functional.output
    output:
        f"{GENESPACE_DIR}/gffs/longest_isoforms.gff3"
    conda: '../envs/annotation.yaml'
    script:
        "../scripts/python/keep_longest_isoform.py"

rule create_subgenome_fasta:
    """
    Creates a fasta with the chromosomes for each subgenome
    """
    input:
        get_create_subgenome_fasta_input
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/split/{ASSEMBLY_NAME}_{{sg}}.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule create_subgenome_gff:
    """
    Create GFF file for each subgenome
    """
    input:
        gff = rules.keep_longest_isoform.output
    output:
        f"{GENESPACE_DIR}/genome_dir/{{sg}}/{{sg}}.gff3"
    shell:
        """
        grep '{wildcards.sg}' {input} > {output}
        """

rule get_subgenome_proteins:
    """
    Extract proteins for each subgenome
    """
    input:
        gff = rules.create_subgenome_gff.output,
        ref = rules.create_subgenome_fasta.output
    output:
        f"{GENESPACE_DIR}/genome_dir/{{sg}}/{{sg}}_proteins.fasta"
    container: 'docker://teambraker/braker3:v.1.0.0'
    shell:
        """
        gffread -y {output} -g {input.ref} {input.gff}
        """

rule genespace:
    input:
        flag = expand(rules.get_subgenome_proteins.output, sg=['Pall', 'Occ']),
        dir = f"{GENESPACE_DIR}/genome_dir"
    output:
        directory(f"{GENESPACE_DIR}/output")
    container: 'library://james-s-santangelo/genespace/genespace:1.2.3'
    threads: 24
    script:
        "../scripts/r/genespace.R"

rule genespace_done:
    input:
        rules.genespace.output
    output:
        f"{GENESPACE_DIR}/genespace.done"
    shell:
        """
        touch {output}
        """
