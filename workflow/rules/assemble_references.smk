# Rules to create haploid reference assembly and phased diploid assembly for annotation and publishing.

rule create_reference_assemblies:
    input:
        rules.fix_haplotypes.output,
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeMapping.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_chromosomesToReverseComplement.csv',
        expand(rules.minimap_hap_vs_hap_paf.output, ver='revised', hap_comp=HAP_VS_HAP)
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_haploid_reference.fasta"
    log: LOG_DIR + '/notebooks/generate_assemblies_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/generate_assemblies.py.ipynb"

rule reference_assemblies_done:
    input:
        rules.create_reference_assemblies.output
    output:
        f'{REFERENCE_ASSEMBLIES_DIR}/reference_assemblies.done'
    shell:
        """
        touch {output}
        """
