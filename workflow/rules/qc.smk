# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule quast_haploid_ref:
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

rule qc_done:
    input:
        rules.quast_haploid_ref.output,
        rules.quast_diploid_ref.output
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
