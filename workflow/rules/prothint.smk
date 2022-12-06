#Prothint

rule prothint_run:
	output:
		gff_prot_basic = f"{ANNOTATION}/prothint/prothint.gff",
		gff_evidence = f"{ANNOTATION}/prothint/evidence.gff",
		gff_prot_ready = f"{ANNOTATION}/prothint/prothint_augutus.gff"
	input:
		plant_db = '../../resources/Viridiplantae_db.fasta',
		masked_genome = f"{ANNOTATION}/repeat_masker/{{hap}}/{{hap}}_softMasked.fasta"
	log: LOG + '/prothint/prot_log/prothint_run.log'
	params:
		outputdir = f"{ANNOTATION}/prothint"
	conda:'../env/prothint.yaml'
	shell:
	"""
		cpanm MCE::Mutex threads YAML Thread::Queue Math::Utils \
		prothint.py {input.masked_genome} \
		{input.plant_db} \
		--workdir {params.outputdir} \
	"""

rule prothint_done:
	output:
		f"{ANNOTATION}/prothint.done"
	input:
		expand(rules.prothint_run.output)
	shell:
	"""
		touch {output}
	"""

rule braker2_run:
	output:
		braker_gtf = f"{ANNOTATION}/braker/braker_dcg/braker.gtf",
		hintsfile = f"{ANNOTATION}/braker/braker_dcg/hintsfile.gff"
	input:
		PROT = {rule.prothint_run.output.gff_prot_ready},
		masked_genome = f"{ANNOTATION}/repeat_masker/{{hap}}/{{hap}}_softMasked.fasta",
		Star_Bam = f"{ANNOTATION_DIR}/star/star_align/Aligned.sortByCoord.out.bam"
	log: LOG + '/braker/braker_log/braker_annotate.log'
	params:
		outputdir = f"{ANNOTATION_DIR}/braker/braker_dcg"
	cores:20
	conda:'../env/braker2.yaml'
	shell:
	"""
		cpanm Hash::Merge List::Util MCE::Mutex Module::Load::Conditional Parallel::ForkManager POSIX Scalar::Util::Numeric YAML File::Spec::FunctionsMath::Utils File::HomeDir
		braker.pl --genome {input.masked_genome} \
		--hints {input.PROT} \
		--bam={input.Star_Bam} \
		-- etpmode \
		--softmasking \
		--cores {cores}
	
	""" 

rule braker2_done:
	input:
	expand(rules.braker2_run.output}
	output:
	f"{ANNOTATION}/braker2.done"
	shell:
	"""
	touch {output}
	"""
