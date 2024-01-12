#######################################################
#	Pipeline Parasito (trouver un meilleur nom)
#
#	Auteur : Julien Robert
#	Mail : julien.robert@aphp.fr
#######################################################

#Anaconda environment : fungi
#Python 3.6.13 (3.6 à cause de snakemake)
#snakemake v3.13.3 (conda)
#picard v2.18.29 (conda)
#gatk4 v4.3.0.0 (conda)
#bwa v0.7.17 (conda)
#samtools v1.6 (conda)
#bcftools v1.9 (conda)
#iqtree v2.0.3 (conda)
#fastqc v0.11.9 (conda)
#multiqc v1.14 (pip install)
#trimmomatic v0.39 (conda)

import random
import glob
import re
import os

#fichier de configuration
#configfile: "config.json"

#Paths
inpath = config["general_path"]["INPUT_PATH"]
outpath = config["general_path"]["OUTPUT_PATH"]
fastapath = config["general_path"]["REF_PATH"]
filterGatkGenotypes = config["general_path"]["filterGatkGenotypes_PATH"]
vcftofasta = config["general_path"]["vcfToFasta_PATH"]

#Etape pré-analyse
#Récupération des noms d'échantillons
sample_id = []
sampleName=glob.glob(inpath+"/*.fastq.gz")
for name in sampleName:
	name = name.replace(inpath+"/", '')
	a = re.split('/', name)
	name = a[0]
	name = name.replace('_1.fastq.gz', '')
	sample_id.append(name)

variants = ' --variant '.join(str(outpath + "/" + elem + "/" + elem + ".g.vcf.gz") for elem in sample_id)

seed = random.randint(10000, 99999)

#Récupération du nom du run
run_name=os.path.basename(inpath)

mate_ids = ["R1","R2"]

#Déclaration des sorties
#trimmomatic = expand((output_path+"/{sample_id}/{sample_id}_R1_paired.fq.gz", output_path+"/{sample_id}/{sample_id}_R2_paired.fq.gz", output_path+"/{sample_id}/{sample_id}_R1_unpaired_trim.fq.gz", output_path+"/{sample_id}/{sample_id}_R2_unpaired_trim.fq.gz"), sample_id = sample_ids),
#fastqc = expand(outpath+"/multiqc/{sample_id}_L001_001_fastqc.html", outpath+"/multiqc/{sample_id}_L001_001_fastqc.zip", sample_id=sample_id)
multiqc = (outpath+"/multiqc/multiqc_report.html")
#ubam = expand(outpath+"/{sample_id}/{sample_id}.unaligned.bam", sample_id=sample_id)
#bwa_mem = expand(outpath+"/{sample_id}/{sample_id}_bwa_mem.bam", sample_id=sample_id)
#markilluminaadapters = expand((outpath+"/{sample_id}/{sample_id}_markilluminaadapters.bam", outpath+"/{sample_id}/{sample_id}_markilluminaadapters_metrics.txt"), sample_id = sample_id)
#addorreplacereadgroups = expand((outpath+"/{sample_id}/{sample_id}.readgroups_unmapped.bam"), sample_id = sample_id)
#samtofastq = expand((outpath+"/{sample_id}/{sample_id}.fastq"), sample_id = sample_id)
#bwa_mem2 = expand((outpath+"/{sample_id}/{sample_id}.sam"), sample_id = sample_id)
#sam2bam = expand((outpath+"/{sample_id}/{sample_id}.bam"), sample_id = sample_id)
#mergebamalignment = expand((outpath+"/{sample_id}/{sample_id}.sorted.bam"), sample_id = sample_id)
#samtools_sort = expand((outpath+"/{ref_id}/{ref_id}.sorted.bam"), ref_id = ref_id)
#markduplicates = expand((outpath+"/{sample_id}/{sample_id}.mark_duplicates.bam", outpath+"/{sample_id}/{sample_id}.mark_duplicates.metrics"), sample_id = sample_id)
#reorderbam = expand((outpath+"/{sample_id}/{sample_id}.reordered.bam", outpath+"/{sample_id}/{sample_id}.reordered.bai"), sample_id =sample_id)
#haplotypecaller = expand((outpath+"/{sample_id}/{sample_id}.g.vcf.gz", outpath+"/{sample_id}/{sample_id}.g.vcf.gz.tbi"), sample_id =sample_id)
#combinegvcfs = (outpath+"/combined_gvcfs.vcf.gz", outpath+"/combined_gvcfs.vcf.gz.tbi")
#genotypegvcfs = (outpath+"/genotyped_gvcfs.vcf.gz", outpath+"/genotyped_gvcfs.vcf.gz.tbi")
#hardfiltration = (outpath+"/indels_filtered.vcf.gz", outpath+"/snps_filtered.vcf.gz", outpath+"/"+run_name+".hard_filtered.vcf.gz)
#makereadablevcf = (outpath+"/"+run_name+".nonbinary.vcf")
#costumevcffilter = (outpath+"/"+run_name+".filtered_SNPs_GQ50_AD08_DP10.vcf", outpath+"/"+run_name+".variant_qc_genotype_filter.tsv")
#vcftofasta = (outpath+"/"+run_name+".vcftofasta.fasta")
#mafft = (outpath+"/"+run_name+".local_alignment.fasta")
iqtree = (outpath+"/"+run_name+".treefile")

rule all:
	input:
#		trimmomatic,
#		fastqc,
		multiqc,
#		ubam,
#		bwa_mem,
#		markilluminaadapters,
#		addorreplacereadgroups,
#		samtofastq,
#		bwa_mem2,
#		sam2bam,
#		mergebamalignment,
#		samtools_sort,
#		markduplicates,
#		reorderbam,
#		haplotypecaller,
#		combinegvcfs,
#		genotypegvcfs,
#		hardfiltration,
#		makereadablevcf,
#		costumevcffilter,
#		vcftofasta,
#		mafft,
		iqtree
	shell:
		"touch "+outpath+"/done"

rule get_sample:
	

rule trimmomatic:
	output: 
		R1 = outpath+"/{sample_id}/{sample_id}_R1_trim.fastq.gz",
		R2 = outpath+"/{sample_id}/{sample_id}_R2_trim.fastq.gz",
		R1u = outpath+"/{sample_id}/{sample_id}_R1_unpaired_trim.fastq.gz",
		R2u = outpath+"/{sample_id}/{sample_id}_R2_unpaired_trim.fastq.gz",
	log:
		outpath+"/logs/{sample_id}_trimmomatic.log"
	params:
		trim_options = config["trimmomatic"]["OPTIONS"],
		trim_params = config["trimmomatic"]["PARAMS"],
	shell:
		'trimmomatic {params.trim_options} {input.R1} {input.R2} {output.R1} {output.R1u} {output.R2} {output.R2u} {params.trim_params}'

rule fastqc:
	input:
		R1 = outpath+"/{sample_id}/{sample_id}_R1_trim.fastq.gz",
		R2 = outpath+"/{sample_id}/{sample_id}_R2_trim.fastq.gz"
	output:
		htmlfileR1 = outpath+"/multiqc/{sample_id}_R1_trim_fastqc.html",
		zipfileR1 = outpath+"/multiqc/{sample_id}_R1_trim_fastqc.zip",
		htmlfileR2 = outpath+"/multiqc/{sample_id}_R2_trim_fastqc.html",
		zipfileR2 = outpath+"/multiqc/{sample_id}_R2_trim_fastqc.zip"
	log : 
		outpath+"/logs/{sample_id}.fastqc.log"
	
	shell :
		'fastqc -o {outpath}/multiqc {input.R1} {input.R2} 2> {log}'

"""
rule trimmomatic:
	output: 
		outpath+"/{sample_id}/{sample_id}_R1_paired.fastq.gz"
	log:
		outpath+"/logs/{sample_id}_trimmomatic.log"
	params:
		trim_options = config["trimmomatic"]["OPTIONS"],
		trim_params = config["trimmomatic"]["PARAMS"],
	shell:
		'trimmomatic {params.trim_options} {inpath}/{wildcards.sample_id}_L001_001.fastq.gz {output} {params.trim_params} 2> {log}'

rule fastqc:
	output:
		html = outpath+"/multiqc/{sample_id}_L001_001_fastqc.html",
		zip = outpath+"/multiqc/{sample_id}_L001_001_fastqc.zip"
	log:
		outpath+"/logs/{sample_id}.fastqc.log"
	
	shell :
		'fastqc -o {outpath}/multiqc {inpath}/{wildcards.sample_id}_L001_001.fastq.gz 2> {log}'
"""
rule multiqc:
	input:
		html = expand(outpath+"/multiqc/{sample_id}_L001_001_fastqc.html", sample_id = sample_id),
		zip = expand(outpath+"/multiqc/{sample_id}_L001_001_fastqc.zip", sample_id = sample_id)
	output:
		outpath+"/multiqc/multiqc_report.html"
	log:
		outpath+"/logs/multiqc.log"
	
	shell :
		'multiqc -f -o {outpath}/multiqc/ {outpath}/multiqc/ --filename multiqc_report.html 2> {log}'

rule ubam:
	input:
		R1 = outpath+"/{sample_id}/{sample_id}_R1_trim.fastq.gz",
		R2 = outpath+"/{sample_id}/{sample_id}_R2_trim.fastq.gz"
	output:
		temp(outpath+"/{sample_id}/{sample_id}.unaligned.bam")
	log:
		outpath+"/logs/{sample_id}.fastqtosam.log"
	shell:
		'picard FastqToSam FASTQ={input.R1} FASTQ2={input.R2} OUTPUT={output} READ_GROUP_NAME=1 SAMPLE_NAME={wildcards.sample_id} LIBRARY_NAME=lib_{wildcards.sample_id} PLATFORM_UNIT=unit1 PLATFORM=ILLUMINA'

"""
rule ubam:
	output:
		temp(outpath+"/{sample_id}/{sample_id}.unaligned.bam")
	log:
		outpath+"/logs/{sample_id}.fastqtosam.log"
	shell:
		'picard FastqToSam FASTQ={inpath}/{wildcards.sample_id}_L001_001.fastq.gz OUTPUT={output} READ_GROUP_NAME=1 SAMPLE_NAME={wildcards.sample_id} LIBRARY_NAME=lib_{wildcards.sample_id} PLATFORM_UNIT=unit1 PLATFORM=ILLUMINA'

rule bwa_mem:
	output:
		temp(outpath+"/{sample_id}/{sample_id}_align.bam")
	params:
		bwamem_db = config["bwamem"]["DB"],
		options = config["bwamem"]["OPTIONS"]
	log:
		outpath+"/logs/{sample_id}.bwamem.log"
	shell:
		'bwa mem {params.options} {params.bwamem_db} {inpath}/{wildcards.sample_id}_L001_001.fastq.gz 2> {log} > {output}'
"""

rule bwa_mem : 
	input : 
		R1 = outpath+"{sample_id}/{sample_id}_R1_trim.fastq.gz",
		R2 = outpath+"{sample_id}/{sample_id}_R2_trim.fastq.gz"
	output : 
		outpath+"{sample_id}/{sample_id}_align.sam"
	params:
		bwamem_db = config["bwamem"]["DB"]
	log :
		"logs/{sample_id}.bwamem.log"
	shell : 
		'bwa mem {params.options} {params.bwamem_db} {input.R1} {input.R2} 2> {log} > {output}'

rule markilluminaadapters:
	input:
		outpath+"/{sample_id}/{sample_id}_align.sam"
	output:
		bam = temp(outpath+"/{sample_id}/{sample_id}_markilluminaadapters.bam"),
		metrics_adapters = outpath+"/{sample_id}/{sample_id}_markilluminaadapters_metrics.txt"
	log:
		outpath+"/logs/{sample_id}.markilluminaadapters.log"
	shell:
		'picard MarkIlluminaAdapters INPUT={input} OUTPUT={output.bam} METRICS={output.metrics_adapters} 2> {log}'

rule addorreplacereadgroups:
	input:
		outpath+"/{sample_id}/{sample_id}.unaligned.bam"
	output:
		temp(outpath+"/{sample_id}/{sample_id}.readgroups_unmapped.bam")
	log:
		outpath+"/logs/{sample_id}.addorreplacereadgroups.log"
	shell: 
		'picard AddOrReplaceReadGroups I={input} O={outpath}/{wildcards.sample_id}/{wildcards.sample_id}.readgroups_unmapped.bam ID=FLOWCELL_{wildcards.sample_id} LB=LIB_{wildcards.sample_id} PL=ILLUMINA SM={wildcards.sample_id} PU=unit1 2> {log}'

rule samtofastq:
	input:
		outpath+"/{sample_id}/{sample_id}_markilluminaadapters.bam"
	output:
		temp(outpath+"/{sample_id}/{sample_id}.fastq")
	log:
		outpath+"/logs/{sample_id}.samtofastq.log"
	shell:
		'picard SamToFastq I={input} FASTQ={output} CLIP_ATTR=XT CLIP_ACT=2 INTER=true NON_PF=true 2> {log}'

rule bwa_mem2:
	input:
		outpath+"/{sample_id}/{sample_id}.fastq"
	output:
		temp(outpath+"/{sample_id}/{sample_id}.sam")
	params:
		ref = config["bwamem"]["DB"],
		options = config["bwamem"]["OPTIONS"]
	log:
		outpath+"/logs/{sample_id}.bwa_mem2.log"
	shell:
		"bwa mem {params.options} -M -R '@RG\\tID:FLOWCELL_{wildcards.sample_id}\\tSM:{wildcards.sample_id}\\tPL:ILLUMINA\\tLB:LIB_{wildcards.sample_id}' -p {params.ref} {input} 2> {log} > {output}"

rule sam2bam:
	input:
		"{prefix}.sam"
	output:
		temp("{prefix}.bam")
	log:
		"{prefix}.sam2bam.log"
	shell:
		'samtools view -1 {input} > {output}'

rule mergebamalignment:
	input:
		ubam = outpath+"/{sample_id}/{sample_id}.readgroups_unmapped.bam",
		bam = outpath+"/{sample_id}/{sample_id}.bam"
	output:
		temp(outpath+"/{sample_id}/{sample_id}.sorted.bam")
	params:
		ref = config["bwamem"]["DB"]
	log:
		outpath+"/logs/{sample_id}.mergebamalignment.log"
	shell:
		'picard MergeBamAlignment ALIGNED={input.bam} UNMAPPED={input.ubam} O={output} R={params.ref} 2> {log}'

rule samtools_sort:
	input:
		outpath+"/{ref_id}/{ref_id}_align_fasta.bam"
	output:
		temp(outpath+"/{ref_id}/{ref_id}.sorted_fasta.bam")
	log:
		outpath+"/logs/{ref_id}.sort.log"
	shell:
		'samtools sort -O bam -o {output} {input} 2> {log}'

rule markduplicates:
	input:
		outpath+"/{sample_id}/{sample_id}.sorted.bam"
	output:
		bam = temp(outpath+"/{sample_id}/{sample_id}.marked_duplicates.bam"),
		metrics_duplicates = outpath+"/{sample_id}/{sample_id}.marked_duplicates.metrics"
	log:
		outpath+"/logs/{sample_id}.markduplicates.log"
	shell:
		'picard MarkDuplicates I={input} O={output.bam} M={output.metrics_duplicates} 2> {log}'

rule reorderbam:
	input:
		outpath+"/{sample_id}/{sample_id}.marked_duplicates.bam"
	output:
		bam = outpath+"/{sample_id}/{sample_id}.reordered.bam",
		bai = outpath+"/{sample_id}/{sample_id}.reordered.bai"
	params:
		ref = config["bwamem"]["DB"]
	log:
		outpath+"/logs/{sample_id}.reorderbam.log"
	run:
		shell('picard ReorderSam INPUT={input} OUTPUT={output.bam} REFERENCE={params.ref} 2> {log}')
		shell('picard BuildBamIndex I={output.bam}')

rule haplotypecaller:
	input:
		bam = outpath+"/{sample_id}/{sample_id}.reordered.bam",
		bai = outpath+"/{sample_id}/{sample_id}.reordered.bai"
	output:
		gvcf = temp(outpath+"/{sample_id}/{sample_id}.g.vcf.gz"),
		gvcf_index = temp(outpath+"/{sample_id}/{sample_id}.g.vcf.gz.tbi")
	params:
		ref = config["bwamem"]["DB"]
	log:
		outpath+"/logs/{sample_id}.haplotypecaller.log"
	shell:
		'gatk HaplotypeCaller -R {params.ref} -I {input.bam} -O {output.gvcf} -ERC GVCF -ploidy 1 2> {log}'


rule combinegvcfs: #rassemblement de tous les échantillons
	input:
		gvcf = expand(outpath+"/{sample_id}/{sample_id}.g.vcf.gz", sample_id = sample_id),
		gvcf_index = expand(outpath+"/{sample_id}/{sample_id}.g.vcf.gz.tbi", sample_id = sample_id)
	output:
		gvcf = temp(outpath+"/combined_gvcfs.vcf.gz"),
		gvcf_index = temp(outpath+"/combined_gvcfs.vcf.gz.tbi")
	params:
		ref = config["bwamem"]["DB"]
	log:
		outpath+"/logs/combinegvcfs.log"
	shell:
		'gatk CombineGVCFs -R {params.ref} -O {output.gvcf} --variant {variants} 2> {log}'

rule genotypegvcfs:
	input:
		gvcf = outpath+"/combined_gvcfs.vcf.gz",
		gvcf_index = outpath+"/combined_gvcfs.vcf.gz.tbi"
	output:
		gvcf = temp(outpath+"/genotyped_gvcfs.vcf.gz"),
		gvcf_index = temp(outpath+"/genotyped_gvcfs.vcf.gz.tbi")
	params:
		ref = config["bwamem"]["DB"]
	log:
		outpath+"/logs/genotypegvcfs.log"
	shell:
		'gatk GenotypeGVCFs -R {params.ref} -O {output.gvcf} -V {input.gvcf} 2> {log}'

rule hardfiltration:
	input:
		gvcf = outpath+"/genotyped_gvcfs.vcf.gz",
		gvcf_index = outpath+"/genotyped_gvcfs.vcf.gz.tbi"
	output:
		indels = temp(outpath+"/indels_filtered.vcf.gz"),
		snps = temp(outpath+"/snps_filtered.vcf.gz"),
		all = temp(outpath+"/"+run_name+".hard_filtered.vcf.gz")
	params:
		ref = config["bwamem"]["DB"]
	log:
		outpath+"/logs/hardfiltration.log"
	run:
		shell('gatk SelectVariants -V {input.gvcf} -R {params.ref} -select-type SNP -O {outpath}/snps.vcf.gz')
		shell('gatk SelectVariants -V {input.gvcf} -R {params.ref} -select-type INDEL -O {outpath}/indels.vcf.gz')
		shell('gatk VariantFiltration -V {outpath}/snps.vcf.gz -R {params.ref} -filter "QD < 20.0" --filter-name "QD20" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -O {output.snps}')
		shell('gatk VariantFiltration -V {outpath}/indels.vcf.gz -R {params.ref} -filter "QD < 20.0" --filter-name "QD20" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -O {output.indels}')
		shell('gatk MergeVcfs -R {params.ref} -I {output.snps} -I {output.indels} -O {output.all}')

rule makereadablevcf:
	input:
		outpath+"/"+run_name+".hard_filtered.vcf.gz"
	output:
		temp(outpath+"/"+run_name+".nonbinary.vcf")
	log:
		outpath+"/logs/makereadablevcf.log"
	shell:
		'bcftools view {input} > {output}'

rule costumevcffilter:
	input:
		outpath+"/"+run_name+".nonbinary.vcf"
	output:
		vcf_filtered = outpath+"/"+run_name+".filtered_SNPs_GQ50_AD08_DP10.vcf",
		vcf_filtered_stats = outpath+"/"+run_name+".variant_qc_genotype_filter.tsv"
	log:
		outpath+"/logs/costumevcffilter.log"
	shell:
		'python {filterGatkGenotypes}/filterGatkGenotypes.py --min_GQ 50 --min_percent_alt_in_AD 0.8 --min_total_DP 10 {input} > {output.vcf_filtered} 2> {output.vcf_filtered_stats}'

rule vcftofasta:
	input:
		outpath+"/"+run_name+".filtered_SNPs_GQ50_AD08_DP10.vcf"
	output:
		outpath+"/"+run_name+".vcftofasta.fasta"
	log:
		outpath+"/logs/vcfToFasta.log"
	run:
		shell('echo {input} > {outpath}/name.txt')
		shell('python {vcftofasta}/vcfSnpsToFasta.py {outpath}/name.txt > {output}')

rule iqtree:
	input:
		outpath+"/"+run_name+".vcftofasta.fasta"
	output:
		outpath+"/"+run_name+".treefile"
	params:
		options = config["iqtree"]["OPTIONS"]
	log:
		outpath+"/logs/iqtree.log"
	run:
		shell("sed 's/\*/N/g' {input} > {outpath}/data_corrected.fasta")
		shell('iqtree -s {outpath}/data_corrected.fasta {params.options} -pre {outpath}/{run_name} 2> {log}')
		#shell('echo {seed} > {outpath}/seed.txt')
