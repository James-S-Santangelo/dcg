#STAR 

rule build_star:
	output:
	 	CHR_LENGTH = f"{ANNOTATION_DIR}/star/star_build/chrLength.txt",
		CHR_NAMELENGTH = f"{ANNOTATION_DIR}/star/star_build/chrNameLength.txt",
		CHR_NAME = f"{ANNOTATION_DIR}/star/star_build/chrName.txt",
		CHR_START = f"{ANNOTATION_DIR}/star/star_build/chrStart.txt",
		GENOME_PARAMS = f"{ANNOTATION_DIR}/star/star_build/genomeParameters.txt",
		SA = f"{ANNOTATION_DIR}/star/star_build/SA",
		SAINDEX = f"{ANNOTATION_DIR}/star/star_build/SAindex",
		GENOME_FASTA = f"{ANNOTATION_DIR}/star/star_build/",
		GENOME_GTF = f"{ANNOTATION_DIR}/star/star_build/",
		WC_GENOME = f"{ANNOTATION_DIR}/star/star_build/Genome"
		
	input:
		masked_genome = f"{ANNOTATION_DIR}/repeat_masker/{{hap}}/{{hap}}_softMasked.fast"
	log: LOG_DIR + '/star/star_log/build.log'
	conda:'../envs/star.yaml'
	threads: 16
	params:
		outdir = f"{ANNOTATION_DIR}/star/star_build"
	shell:
		"""
		STAR --runMode genomeGenerate \
		--genomeDir {params.outdir} \
		--genomeFastaFile {input.masked_genome} \
		--runThreadN {threads}
		"""

rule build_done:
	input:
		expand(rules.build_star.output)
	output:
		f"{ANNOTATION_DIR}/build_star.done"
	shell:
		"""
		touch {output}
		"""
rule align_star:
	output:
		star_align = f"{ANNOTATION_DIR}/star/star_align/Aligned.sortedByCoord.out.bam"	
	input:
		R1fq = f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fq",
		R2fq = f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fq",
		STAR_BUILD = {rules.params.outdir}
		 
	log: LOG + '/star/star_log/star_align.log'
	conda:'../envs/star.yaml'
	shell:
		"""
		STAR --readFilesIn {input.R1fq} {input.R2fq} \
		--alignIntronMin 20 \
		--alignIntronMax 500000 \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
		--genomeDir {input.STAR_BUILD}  
		"""

rule align_done:
	input:
		expand(rules.align_star.output)
	output:
		f"{ANNOTATION_DIR}/align.done"
	shell:
		"""
		touch {output}
		"""
