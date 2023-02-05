# Snakefile with rules to fix haplotype issues identified through manual inspection of dotplots

rule fix_haplotypes:
    input:
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.agp" for hap in HAPS],
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.fasta" for hap in HAPS],
        expand(rules.dotplot_hap_vs_hap.output, ver='original', hap_comp=HAP_VS_HAP)
    output:
        expand(f"{REVISED_HAP_DIR}/{{hap}}_revised.fasta", hap=HAPS)
    log: LOG_DIR + '/notebooks/fix_haplotypes_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/fix_haplotypes.py.ipynb"

rule revised_haplotypes_done:
    input:
        rules.fix_haplotypes.output
    output:
        f'{REVISED_HAP_DIR}/revised_haplotypes.done'
    shell:
        """
        touch {output}
        """
