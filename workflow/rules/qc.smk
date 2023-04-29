# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule quast_haploid_ref:
    """
    Use QUAST to QC the haploid reference assembly
    """
    input:
        fasta = rules.repeat_masker.output.fasta,
        gff = f"{rules.funannotate_annotate.output}/annotate_results/Trifolium_repens.gff3"
    output:
        directory(f"{QC_DIR}/quast/haploid")
    log: LOG_DIR + '/quast/quast_haploid.log'
    conda: '../envs/qc.yaml'
    threads: 8
    shell:
        """
        quast.py -o {output} \
            --features {input.gff} \
            --threads {threads} \
            --split-scaffolds \
            --large \
            --k-mer-stats \
            {input.fasta} &> {log}
        """

rule quast_diploid_ref:
    """
    Use QUAST do QC the phased diploid assembly.
    """
    input:
        rules.create_reference_assemblies.output[1],
        rules.create_reference_assemblies.output[2]
    output:
        directory(f"{QC_DIR}/quast/diploid")
    log: LOG_DIR + '/quast/quast_diploid.log'
    conda: '../envs/qc.yaml'
    threads: 8
    shell:
        """
        quast.py -o {output} \
            --threads {threads} \
            --split-scaffolds \
            --large \
            --k-mer-stats \
            {input} &> {log}
        """

rule run_busco:
    """
    Assess annotation completeness by running BUSCO against both the embryophyta and Fabales databases.
    """
    input:
        prot = f"{rules.funannotate_annotate.output}/annotate_results/Trifolium_repens.proteins.fa"
    output:
        directory(f"{QC_DIR}/busco/TrR_v6_{{db}}")
    log: LOG_DIR + '/busco/busco_{db}.log'
    conda: '../envs/qc.yaml'
    threads: 32
    params:
        out_path = f"{QC_DIR}/busco/",
        out_name = "TrR_v6_{db}"
    shell:
        """
        busco -m protein \
            -i {input.prot} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            --lineage {wildcards.db} \
            --force \
            --cpu {threads} &> {log}
        """

rule qc_done:
    input:
        rules.quast_haploid_ref.output,
        rules.quast_diploid_ref.output,
        expand(rules.run_busco.output, db = ['embryophyta_odb10', 'fabales_odb10'])
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
