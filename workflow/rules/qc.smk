# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule quast_haploid_ref:
    """
    Use QUAST to QC the haploid reference assembly
    """
    input:
        fasta = rules.repeat_masker.output.fasta,
        gff = rules.funannotate_annotate.output.gff3
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
        rules.create_reference_assemblies.output.dip_hap1,
        rules.create_reference_assemblies.output.dip_hap2
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

rule run_busco_protein:
    """
    Assess annotation completeness by running BUSCO against both the embryophyta and Fabales databases.
    """
    input:
        prot = rules.get_final_proteins.output 
    output:
        directory(f"{QC_DIR}/busco/{ASSEMBLY_NAME}_{{db}}")
    log: LOG_DIR + '/busco/busco_{db}.log'
    conda: '../envs/qc.yaml'
    threads: 24
    params:
        out_path = f"{QC_DIR}/busco/",
        out_name = f"{ASSEMBLY_NAME}_{{db}}"
    shell:
        """
        busco -m protein \
            -i {input} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            --lineage {wildcards.db} \
            --force \
            --cpu {threads} &> {log}
        """

rule run_busco_genome:
    input:
        lambda w: rules.repeat_masker.output.fasta if w.ass == 'hap' else (rules.create_reference_assemblies.output.dip_hap1 if w.ass == 'dip1' else rules.create_reference_assemblies.output.dip_hap2)
    output:
        directory(f"{QC_DIR}/busco/{ASSEMBLY_NAME}_{{ass}}_{{db}}")
    log: LOG_DIR + '/busco/busco_{{ass}}_{db}.log'
    conda: '../envs/qc.yaml'
    threads: 16
    params:
        out_path = f"{QC_DIR}/busco/",
        out_name = f"{ASSEMBLY_NAME}_{{ass}}_{{db}}"
    shell:
        """
        busco -m genome\
            -i {input} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            --lineage {wildcards.db} \
            --force \
            --cpu {threads} &> {log}
        """

rule functional_stats:
    """
    Use AGAT to generate summary statistics of functional annotation
    """
    input:
        gff = rules.gff_sort_functional.output,
        fasta = rules.repeat_masker.output.fasta
    output:
        directory(f"{QC_DIR}/agat_funcStats")
    log: f"{LOG_DIR}/agat_funcStats/agat_funcStats.log"
    container: 'docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0'
    shell:
        """
        agat_sp_functional_statistics.pl --gff {input.gff} \
            -gs {input.fasta} \
            --output {output} 2> {log}
        """

rule qc_done:
    input:
        rules.functional_stats.output,
        rules.quast_haploid_ref.output,
        rules.quast_diploid_ref.output,
        expand(rules.run_busco_protein.output, db = ['embryophyta_odb10', 'fabales_odb10']),
        expand(rules.run_busco_genome.output, ass = ['hap', 'dip1', 'dip2'], db = ['embryophyta_odb10', 'fabales_odb10'])
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
