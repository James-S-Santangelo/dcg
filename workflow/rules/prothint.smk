#Prothint

rule viridiplantae_orthodb:
	output:
		Plant_ProrteinDB = "f{ANNOTATION_DIR}/orthodb/Viridiplantae_protein.fasta"
	log: LOG_DIR + "f{ANNOTATION_DIR}/orthodb/ortho.log"
	params:
		outdir = temp"f{ANNOTATION_DIR}/orthodb"
	shell:
	"""
		wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz \
		tar -zxf odb10_plants_fasta.tar.gz -C {params.outdir}
		cat {params.outdir}/plants/Rawdata/*" > {output}
	"""

rule prothint_run:
	output:
		directory(f"{ANNOTATION_DIR}/prothint"
	input:
		plant_db = rule.viridiplantae_orthodb.Plant_ProteinDB,
		masked_genome = f"{ANNOTATION}/repeat_masker/{{hap}}/{{hap}}_softMasked.fasta"
	log: LOG + '/prothint/prot_log/prothint_run.log'
	params:
		outputdir = f"{ANNOTATION}/prothint"
	container:'singularity://james-s-santangelo/braker/braker'
	shell:
	"""
		prothint.py {input.masked_genome} \
		{input.plant_db} \
		--workdir {params.outputdir} \
	"""

rule braker_rnaseq:
	output:
		hints_rna = f{"ANNOTATION}/braker/braker1_out/hintsfile.gff"	
		aug_hint_rna = f"{ANNOTATION}/braker/braker1_out/augustus.hints.gtf"
	input:
		masked_genome = f"{ANNOTATION}/repeat_masker/{{hap}}/{{hap}}_softMasked.fasta",
		Star_Bam = f"{ANNOTATION_DIR}/star/star_align/{{hap}}_{{acc}}_sorted.bam"
	log: LOG_DIR + '/braker/braker_log/braker_annotate.log'
	params:
		outputdir = f"{ANNOTATION_DIR}/braker/braker1_out"
	cores:20
	conda:'singularity://james-s-santangelo/braker/braker'
	shell:
	"""
		braker.pl --genome {input.masked_genome} \
		--bam={input.Star_Bam} \
		--softmasking \
		--cores {cores} \
		--workingdir={params.outputdir} \
		--species="Trifolium repens" \
		--AUGUSTUS_hints_preds=augustus.hints.gtf
	""" 

rule braker_protein:
	output:
		hints_protein = f{"ANNOTATION}/braker/braker2_out/hintsfile.gff"
		aug_hint_protein = f"{ANNOTATION}/braker/braker2_out/augustus.hints.gtf"
	input:
		PROT = "f{ANNOTATION_DIR}/prothint/prothint_augustus.gff",
		masked_genome = f"{ANNOTATION}/repeat_masker/{{hap}}/{{hap}}_softMasked.fasta",
	log: LOG_DIR + '/braker/braker_log/braker_annotate.log'
	params:
		outputdir = f"{ANNOTATION_DIR}/braker/braker2_out"
	cores:20
	conda:'singularity://james-s-santangelo/braker/braker'
	shell:
	"""
		braker.pl --genome {input.masked_genome} \
		--hints {input.PROT} \
		--epmode \
		--softmasking \
		--cores {cores} \
		--workingdir={params.outputdir} \
		--species="Trifolium repens" \
		--AUGUSTUS_hints_preds=augustus.hints.gtf
	""" 

#COMBINE BRAKER RNA-SEQ with BRAKER PROTEIN

rule tsebra_combine:
	input:
		rna_aug = {rule.braker_rnaseq.output.aug_hint_rna}
		protein_aug = {rule.braker_protein.output.aug_hint_protein} 
		hints_rna = {rule.braker_rnaseq.output.hints_rna}
		hints_protein = {rule.braker_protein.output.hints_protein}
	output:
		braker_combined = f"{ANNOTATION_DIR}/braker/tsebra_combined/braker1+2_combined.gft"
	log: LOG_DIR + '/braker/braker_log/tsebra.log'
	params:
		outdir = f"{ANNOTATION_DIR}/braker/tsebra_combined"
	shell:
	"""
		./bin/tesbra.py -g {input.rna_aug},{input.protein_aug} \
		-c default.cfg \
		-e {input.hints_rna},{input.hints_protein} \ 
		-o f"{params.outdir}/braker1+2_combined.gft"
	"""	
