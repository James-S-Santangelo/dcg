# Rules to create haploid reference assembly and phased diploid assembly for annotation and publishing.

rule fix_haplotypes:
    """
    Imposes manual fixes of raw Dovetail haplotypes and generates FASTA files with revised haplotype sequences. 
    """
    input:
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.agp" for hap in HAPS],
        [HAPLOTYPE_FASTA_DIR + f"/{hap}.fasta" for hap in HAPS],
        expand(rules.dotplot_hap_vs_hap.output, ver='original', hap_comp=HAP_VS_HAP)
    output:
        expand(f"{REVISED_HAP_DIR}/{{hap}}_revised.fasta", hap=HAPS)
    log: LOG_DIR + '/notebooks/fix_haplotypes_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/fix_assembly_issues.py.ipynb"

rule identify_chromosomes_linkageGroups:
    """
    Determines which of the revised haplotypes correspond to which chromosome based on (1) minimap alignments to previous reference genome from Griffiths et al. (2019) and (2) blast alignments of linkage markers from Olsen et. al (2022) F2 mapping population to each haplotype. 

    Produces a number of outputs used to later generate the reference assemblies including: (1) Plot of number of BLAST hits for each linkage group and subgenome, (2) CSV file with lengths of each chromosome in new haplotypes, (3) CSV file with longest chromosome from each of the two haplotypes for each chromosome, (4) CSV file with mappings of chromosomes in the new haploid assembly to those in the previous assembly, (5) Plot of alignments to previous assembly to determine correct orientation of chromosomes, (6) CSV file with list of chromosomes that need to be reverse complemented to match previous assembly. 
        """
    input:
        expand(rules.blast_linkageMarkers.output, hap=HAPS, map_pop=['SG','DG']),
        f"{GENMAP_RESOURCE_DIR}/DG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_marker_key.csv",
        f"{GENMAP_RESOURCE_DIR}/DG_genMap.csv",
        f"{GENMAP_RESOURCE_DIR}/SG_genMap.csv",
        TRR_FIVE_SCAFFS,
        expand(rules.minimap_hap_vs_TrRvFive.output, hap=HAPS)
    output:
        blast = f'{REFERENCE_ASSEMBLIES_DIR}/{ASSEMBLY_NAME}_vs_TrR5_and_LG/BLASThits_linkageGroup_by_subGenome.pdf',
        scaf_len = f'{REFERENCE_ASSEMBLIES_DIR}/{ASSEMBLY_NAME}_vs_TrR5_and_LG/{ASSEMBLY_NAME}_scaffoldLengths_bySubGenome_and_Hap.csv',
        logScaf_LG = f'{REFERENCE_ASSEMBLIES_DIR}/{ASSEMBLY_NAME}_vs_TrR5_and_LG/{ASSEMBLY_NAME}_longestScaffolds_byLinkageGroup.csv',
        chrom_map = f'{REFERENCE_ASSEMBLIES_DIR}/{ASSEMBLY_NAME}_vs_TrR5_and_LG/TrR_v5_to_{ASSEMBLY_NAME}_chromosomeMapping.csv',
        chrom_or = f'{REFERENCE_ASSEMBLIES_DIR}/{ASSEMBLY_NAME}_vs_TrR5_and_LG/TrR_v5_to_{ASSEMBLY_NAME}_chromosomeOrientation.pdf',
        to_revComp = f'{REFERENCE_ASSEMBLIES_DIR}/{ASSEMBLY_NAME}_vs_TrR5_and_LG/{ASSEMBLY_NAME}_chromosomesToReverseComplement.csv',
    log: LOG_DIR + '/notebooks/TrRvSix_vs_TrRvFive_and_LG_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/identify_chromosomes_linkageGroups.r.ipynb"

rule create_reference_assemblies:
    """
    Creates haploid mapping reference assembly and phased diploid assembly from revised haplotype sequences.
    """
    input:
        rules.fix_haplotypes.output,
        rules.identify_chromosomes_linkageGroups.output.chrom_map,
        rules.identify_chromosomes_linkageGroups.output.to_revComp,
        rules.identify_chromosomes_linkageGroups.output.scaf_len
    output:
        haploid = f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/{ASSEMBLY_NAME}_haploid_reference.fasta",
        dip_hap1 = f"{REFERENCE_ASSEMBLIES_DIR}/diploid_reference/{ASSEMBLY_NAME}_hap1_reference.fasta",
        dip_hap2 = f"{REFERENCE_ASSEMBLIES_DIR}/diploid_reference/{ASSEMBLY_NAME}_hap2_reference.fasta",
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
