configfile: 'config.yaml'

#module load python/3.5.1 samtools freebayes vcftools bwa

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}


def return_filename_by_sampname(sampname):
	filenames = []
	if sample_by_name[sampname]['paired']:
		filenames.append(sample_by_name[sampname]['readsfile1'])
		filenames.append(sample_by_name[sampname]['readsfile2'])
	else:
		filenames.append(sample_by_name[sampname]['readsfile'])
	return filenames

def return_file_relpath_by_sampname(wildcards):
	sampname = wildcards.samplename
	pathprefix = sample_by_name[sampname]["path"]
	filesin = return_filename_by_sampname(sampname)
	pathsout = ["".join([pathprefix, fq]) for fq in filesin]
	return pathsout

rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R0.fastq"],
		jason = "meta/FASTP/{samplename}.json"
	params:
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json meta/FASTP/{samplename}.json --html meta/FASTP/{samplename}.html", 
		se_params = "",
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]}"


rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
		jason = "meta/FASTP/{samplename}.json"
	params:
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json meta/FASTP/{samplename}.json --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]}"

#request FASTP reportBacks instead, process them into a unified SQL upload?
# rule clean_all_samps:
# 	input: 
# 		clean_se = [sample_by_name[nom]['path']+nom+".clean.R0.fastq" for nom in sample_by_name.keys() if not sample_by_name[nom]['paired']],
# 		clean_pe = [sample_by_name[nom]['path']+nom+".clean.R"+arr+".fastq" for nom in sample_by_name.keys() if sample_by_name[nom]['paired'] for arr in ["1","2"]],
# 	output:
# 		outflg = "allclean.flag"
# 	shell:
# 		"touch {output.outflg}"

rule FASTP_summarizer:
	input: 
		jason = "meta/FASTP/{samplename}.json"
	output:
		jason_pruned = "meta/FASTP/{samplename}.json.pruned"
	shell:
		"python scripts/fastp_reporter.py {input.jason} {output.jason_pruned} --t {wildcards.samplename}"

rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("meta/FASTP/{samplename}.json.pruned", samplename = sample_by_name.keys())
	output:
		summary = "meta/sequenced_reads.dat"
	run:
		"cat {input.jasons_in} > {output.summary}"


#	request stats/idxstats/flagstats?  
rule bwa_align:
	input:
		reads_in = lambda wildcards: expand("{path}{sample}.clean.R{arr}.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample, arr=[ [1,2] if sample_by_name[wildcards.sample]['paired'] else [0] ][0]),
		ref_genome_file = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]['path'],
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam",
	run:
		shell("bwa aln {input.ref_genome_file} {input.reads_in[0]} > {input.reads_in[0]}.sai ")
		if sample_by_name[wildcards.sample]['paired']:
			shell("bwa aln {input.ref_genome_file} {input.reads_in[1]} > {input.reads_in[1]}.sai ")
			shell("bwa sampe {input.ref_genome_file} {input.reads_in[0]}.sai {input.reads_in[1]}.sai {input.reads_in[0]}  {input.reads_in[1]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		else:
			shell("bwa samse {input.ref_genome_file} {input.reads_in[0]}.sai {input.reads_in[0]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		shell("samtools index {output.bam_out}")


#	request stats/idxstats/flagstats?  

rule bwa_uniq:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam"
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwaUniq.sort.bam"
	params:
		quality="-q 20 -F 0x0100 -F 0x0200 -F 0x0300 -F 0x04",
		uniqueness="XT:A:U.*X0:i:1.*X1:i:0",
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		#original; no dedupe
		#"samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome} - | samtools sort -o {output.bam_out} - "
		shell("samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome_file} - | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - {output.bam_out}")
		shell("samtools index {output.bam_out}")


rule joint_vcf_caller:
	input:
		bams_in = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam", sample=sample_by_name.keys(), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
	output:
		vcf_out = "all_samples.vs_{ref_genome}.{aligner}.vcf"
	params:
		freebayes="--standard-filters",
		runmem_gb=8,
		runtime="12:00:00"
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		shell("freebayes {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.joint_vcf}")

rule write_report:
	input:
		sequenced_reads_summary=["meta/sequenced_reads.dat"],
	output:
		pdf_out="thingy.pdf"
	run:
		"touch {output.pdf_out}"
