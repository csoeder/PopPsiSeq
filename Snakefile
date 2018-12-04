configfile: 'config.yaml'

#module load python/3.5.1 samtools freebayes vcftools bwa bedtools r/3.5.0 rstudio/1.1.453

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


rule reference_genome_reporter:
	input:
		fai_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['fai'],
	output:
		report_out = "meta/reference_genomes/{ref_gen}.fai.report"
	shell:
		"""
		mkdir -p meta/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule demand_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("meta/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "meta/reference_genomes.summary"
	shell:
		"cat {input.refgen_reports} > {output.refgen_summary}"



rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R0.fastq"],
		jason = "{pathprefix}/{samplename}.False.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.False.json",# --html meta/FASTP/{samplename}.html", 
		se_params = "",
	message:
		"FASTP QA/QC on single-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]}"


rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
		jason = "{pathprefix}/{samplename}.True.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.True.json",# --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	message:
		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
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

lambda wildcards: return_file_relpath_by_sampname(wildcards)

rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("{path}{samp}.{pairt}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, pairt = sample_by_name[wildcards.samplename]['paired'])
	output:
		jason_pruned = "meta/FASTP/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		cp {input.jason} meta/FASTP/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""

rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("meta/FASTP/{samplename}.json.pruned", samplename = sample_by_name.keys())
	output:
		summary = "meta/sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"


#	request stats/idxstats/flagstats?  
rule bwa_align:
	input:
		reads_in = lambda wildcards: expand("{path}{sample}.clean.R{arr}.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample, arr=[ [1,2] if sample_by_name[wildcards.sample]['paired'] else [0] ][0]),
		ref_genome_file = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]['path'],
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam",
	params:
		runmem_gb=64,
		runtime="12:00:00",
	message:
		"aligning reads from {wildcards.sample} to reference_genome {wildcards.ref_genome} .... "
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
		runmem_gb=16,
		runtime="6:00:00",
	message:
		"filtering alignment of {wildcards.sample} to {wildcards.ref_genome} for quality and mapping uniqueness.... "	
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		#original; no dedupe
		#"samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome} - | samtools sort -o {output.bam_out} - "
		shell("samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome_file} - | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - {output.bam_out}")
		shell("samtools index {output.bam_out}")



rule bam_reporter:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam"
	output:
		report_out = "meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=8,
		runtime="4:00:00",
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	run:
		ref_genome_idx=ref_genome_by_name[wildcards.ref_genome]['fai']
		shell("samtools idxstats {input.bam_in} > {input.bam_in}.idxstats")
		shell("samtools flagstat {input.bam_in} > {input.bam_in}.flagstat")
		shell("bedtools genomecov -max 1 -ibam {input.bam_in} -g {ref_genome_idx} > {input.bam_in}.genomcov")
		#change the -max flag as needed to set 
		shell("""samtools depth -a {input.bam_in} | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' > {input.bam_in}.dpthStats""")
		#https://www.biostars.org/p/5165/
		#save the depth file and offload the statistics to the bam_summarizer script?? 
		shell("python3 scripts/bam_summarizer.py -f {input.bam_in}.flagstat -i {input.bam_in}.idxstats -g {input.bam_in}.genomcov -d {input.bam_in}.dpthStats -o {output.report_out} -t {wildcards.sample}")


rule demand_BAM_analytics:
	input:
		bam_reports = lambda wildcards: expand("meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sample_by_name.keys(), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
	output:
		full_report = "meta/alignments.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
	message:
		"collecting all alignment metadata.... "
	shell:
		"cat {input.bam_reports} > {output.full_report}"
#	cat test.samtools.idxstats | sed \$d | awk '{print $1, $3/$2}' > per_contig_coverage_depth
#	echo  $( cat test.samtools.idxstats |  sed \$d | cut -f 2 | paste -sd+ | bc) $(cat test.samtools.idxstats |  sed \$d | cut -f 3 | paste -sd+ | bc) | tr  " " "\t" | awk '{print "total\t"$1/$2}'


rule joint_vcf_caller:
	input:
		bams_in = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam", sample=sample_by_name.keys(), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
	output:
		vcf_out = "variants/all_samples.vs_{ref_genome}.{aligner}.vcf"
	params:
		freebayes="--standard-filters",
		runmem_gb=8,
		runtime="12:00:00",
	message:
		"Jointly calling variants from all samples mapped to {wildcards.ref_genome} with {wildcards.aligner}"
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		shell("freebayes {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}")


## VCF summary report? eg, shared fraction of genome covered by x% of samples at mindepth?



rule vcf_reporter:
	input:
		vcf_in = "variants/{prefix}.vs_{ref_genome}.{aligner}.vcf"
	output:
		report_out = "meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=8,
		runtime="4:00:00",
#	message:
#		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	shell:
		"""
		cat {input.vcf_in}  | grep -v "#" | cut -f 1 | sort | uniq -c > {output.report_out}.snpsPerContig.tmp
		cat {output.report_out}.snpsPerContig.tmp | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' | cat {output.report_out}.snpsPerContig.tmp - > {output.report_out}.snpsPerContig
		rm {output.report_out}.snpsPerContig.tmp
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --freq 
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --counts
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-indv
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-site
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --singletons

		ref_genome={wildcards.ref_genome}

		tail -n 1 {output.report_out}.snpsPerContig | awk '{print "total_snp_count\t"$1}' | sed -e 's/^/'$ref_genome'\t/g' > {output.report_out}
		"""
	#cat  all_samples.vs_droSec1.bwaUniq.summary.frq.count| cut -f 3 | tail -n +2 | sort | uniq -c
	#####	bi, tri, and quadralelic counts ^^ 
	#replace some of this with vcftools::vcf-stats ?

rule summon_VCF_analytics_base:
	input:
		bam_reports = lambda wildcards: expand("meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary", prefix=wildcards.prefix, ref_genome=["droSim1","droSec1"], aligner="bwaUniq")
	output:
		full_report = "meta/{prefix}.calledVariants.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
	message:
		"collecting all alignment metadata.... "
	shell:
		"""
		prefix={wildcards.prefix}
		cat {input.bam_reports} | sed -e 's/^/'$prefix'\t/g'> {output.full_report}
		"""

## benchmarking??
#	https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#benchmark-rules



rule write_report:
	input:
		reference_genome_summary = ["meta/reference_genomes.summary"],
		sequenced_reads_summary=["meta/sequenced_reads.dat"],
		alignment_summaries = expand("meta/alignments.vs_{ref_genome}.{aligner}.summary", ref_genome=['droSim1', 'droSec1'], aligner=['bwa','bwaUniq']),
	output:
		pdf_out="thingy.pdf"
	params:
		runmem_gb=8,
		runtime="1:00:00",
	message:
		"writing up the results.... "
	run:
		pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
		pwd = subprocess.check_output("pwd",shell=True).decode().rstrip()+"/"
		shell(""" R -e "setwd('{pwd}');Sys.setenv(RSTUDIO_PANDOC='{pandoc_path}')" -e  "peaDubDee='{pwd}'; rmarkdown::render('scripts/PopPsiSeq_summary.Rmd',output_file='{pwd}{output.pdf_out}')"  """)
#		shell(""" R -e "setwd('{pwd}/');" -e Sys.setenv"(RSTUDIO_PANDOC='{pandoc_path}')" -e  rmarkdown::render"('scripts/PopPsiSeq_summary.Rmd',output_file='{output.pdf_out}')"  """)
#		shell(""" R -e Sys.setenv"(RSTUDIO_PANDOC='{pandoc_path}')" -e  "thetitle='My title'; theauthor='me'; rmarkdown::render('scripts/PopPsiSeq_summary.Rmd',output_file='{output.pdf_out}')"  """)
#		shell(""" cp scripts/{output.pdf_out} {output.pdf_out} """)
#pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
#R -e Sys.setenv"(RSTUDIO_PANDOC='$pandoc_path')" -e  rmarkdown::render"('$markDown_in',output_file='$pdf_Out')"
#R -e "setwd('/proj/cdjones_lab/csoeder/PopPsiSeq')" -e Sys.setenv"(RSTUDIO_PANDOC='$pandoc_path')" -e  rmarkdown::render"('scripts/PopPsiSeq_summary.Rmd',output_file='test.pdf')"






