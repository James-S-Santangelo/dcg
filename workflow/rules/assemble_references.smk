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

    Produces a number of outputs used to later generate the reference assemblies including: (1) Plot of number of BLAST hits for each linkage group and subgenome, (2) CSV file with lengths of each chromosome in new haplotypes, (3) CSV file with longest chromosome from each of the two haplotypes for each chromosome, (4) Plot of alignments to previous assembly to determine correct orientation of chromosomes, (5) CSV file with list of chromosomes that need to be reverse complemented to match previous assembly. 
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
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/BLASThits_linkageGroup_by_subGenome.pdf',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_scaffoldLengths_bySubGenome_and_Hap.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrRv6_longestScaffolds_byLinkageGroup.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeMapping.csv',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v5_to_v6_chromosomeOrientation.pdf',
        f'{QC_DIR}/TrRv6_vs_TrR5_and_LG/TrR_v6_chromosomesToReverseComplement.csv',
    log: LOG_DIR + '/notebooks/TrRvSix_vs_TrRvFive_and_LG_processed.ipynb'
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/TrRvSix_vs_TrRvFive_and_LG.r.ipynb"

rule create_reference_assemblies:
    """
    Creates haploid mapping reference assembly and phased diploid assembly from revised haplotype sequences.
    """
    input:
        rules.fix_haplotypes.output,
        rules.identify_chromosomes_linkageGroups.output[3],
        rules.identify_chromosomes_linkageGroups.output[5],
        rules.identify_chromosomes_linkageGroups.output[1]
    output:
        f"{REFERENCE_ASSEMBLIES_DIR}/haploid_reference/TrR_v6_haploid_reference.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/diploid_reference/TrR_v6_hap1_reference.fasta",
        f"{REFERENCE_ASSEMBLIES_DIR}/diploid_reference/TrR_v6_hap2_reference.fasta",
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
