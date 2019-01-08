configfile: 'config.yaml'

#module load python/3.5.1 samtools freebayes vcftools bwa bedtools r/3.5.0 rstudio/1.1.453

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
sampname_by_group = {}
for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]




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
	params:
		runmem_gb=1,
		runtime="5:00",
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
	params:
		runmem_gb=1,
		runtime="5:00",
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
		runtime="64:00:00",
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
		runmem_gb=32,
		runtime="96:00:00",
	message:
		"Jointly calling variants from all samples mapped to \ {wildcards.ref_genome} \ with \ {wildcards.aligner} \ "
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		shell("freebayes {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}")


## VCF summary report? eg, shared fraction of genome covered by x% of samples at mindepth?



rule vcf_reporter:
	input:
		vcf_in = "variants/{prefix}.vs_{ref_genome}.{aligner}.vcf"
	output:
		report_out = "meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary",
		frq_out = "meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary.frq"
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

		tail -n 1 {output.report_out}.snpsPerContig | awk '{{print "total_snp_count\t"$1}}' | sed -e 's/^/'$ref_genome'\t/g' > {output.report_out}
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


rule subset_VCF_to_subgroup:
	input:
		vcf_in = "variants/{prefix}.vs_{ref_genome}.{aligner}.vcf"
	output:
		vcf_out = "variants/{prefix}.subset_{subgroup}.vs_{ref_genome}.{aligner}.vcf"
	params:
		runmem_gb=8,
		runtime = "1:00:00",
	message:
		"Subsetting the variant file \ {input.vcf_in} \ to the individuals in the subgroup \ {wildcards.subgroup} \ .... "
	run:
		member_list = "%s,"*len(sampname_by_group[wildcards.subgroup]) % tuple(sampname_by_group[wildcards.subgroup])
		shell("vcf-subset {input.vcf_in} -u -c {member_list} | vcftools --min-alleles 2 --max-alleles 2  --max-missing-count 1  --vcf - --recode --stdout > {output.vcf_out}")



rule window_maker:
	output:
		windowed='utils/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed'
	params:
		runmem_gb=8,
		runtime="5:00"
	run:
		fai_path = ref_genome_by_name[wildcards.ref_genome]['fai'],
		shell("mkdir -p utils")
		shell(
			'bedtools makewindows -w {wildcards.window_size} -s {wildcards.slide_rate} -g {fai_path} -i winnum | bedtools sort -i - > {output.windowed}'
		)



#######		windowed frequency comparison between groups 	#######################

rule clean_groupFreqs:
	input:
		#vcf_in = "variants/{prefix}.subset_{grup}.vs_{ref_genome}.{aligner}.vcf",
		frq_in = "meta/VCFs/{prefix}.subset_{grup}.vs_{ref_genome}.{aligner}.summary.frq"
	output:
		frq_out = "variant_analysis/freqs/{prefix}.subset_{grup}.vs_{ref_genome}.{aligner}.frq",
	params:
		runmem_gb=16,
		runtime="15:00"
	shell:
		"""
		mkdir -p variant_analysis/freqs/
		tail -n +2 {input.frq_in}  | awk '{{print $1,$2,$2+1,$4,$5,$6}}' | tr " " "\t" > {output.frq_out}
		"""


rule calc_frq_shift:
	input:
		par1_frq = "variant_analysis/freqs/{prefix}.subset_{grup_par1}.vs_{ref_genome}.{aligner}.frq",
		par2_frq = "variant_analysis/freqs/{prefix}.subset_{grup_par2}.vs_{ref_genome}.{aligner}.frq",
		off_frq = "variant_analysis/freqs/{prefix}.subset_{grup_off}.vs_{ref_genome}.{aligner}.frq",
	output:
		frqShft_out = "variant_analysis/freqShift/{prefix}.{grup_off}_with_{grup_par1}_and_{grup_par2}.vs_{ref_genome}.{aligner}.frqShift"
	params:
		runmem_gb=16,
		runtime="1:00:00"
	shell:
		"""
		mkdir -p variant_analysis/freqShift/
		bedtools intersect -wa -wb -a {input.par1_frq} -b  {input.par2_frq} | bedtools intersect -wa -wb -a - -b  {input.off_frq} | cut -f 1,2,4-6,10,12,16,18 | tr ":" "\t" | awk '{{print $1,$2,$2+1,"0","0","+",$4,$6,$3,$7,$8,$10,$11,$13}}' | tr " " "\t" > {output.frqShft_out}.pre
		Rscript scripts/freqShifter.R {output.frqShft_out}.pre {output.frqShft_out}
		"""
		#Rscript scripts/freqShifter.R {output.frqShft_out}.pre {output.frqShft_out}

#	no need for biallelic check - that's done in the subsetting

#	add a check to make sure that it's the same allele in each case, write down the cases which aren't????



rule window_frq_shift:
	input:
		frqShft_in = "variant_analysis/freqShift/{frqshft_prefix}.vs_{ref_genome}.{aligner}.frqShift",
		windows_in = "utils/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed"
	output:
		windowed_out = "variant_analysis/freqShift/{frqshft_prefix}.vs_{ref_genome}.{aligner}.windowed_w{window_size}_s{slide_rate}.frqShift"
	params:
		runmem_gb=16,
		runtime="1:00:00"
	shell:






#cat utils/droSim1_w100000_s100000.windows.bed | grep -w "chr2L" > dev/droSim1_w100000_s100000.windows.chr2L.bed  

#vcftools  --vcf PopSech.chr2L.vs_droSim1.bwaUniq.vcf --stdout --freq 
#vcftools  --vcf PopSim.chr2L.vs_droSim1.bwaUniq.vcf --out PopSim.chr2L.vs_droSim1.bwaUniq --freq 
#vcftools  --vcf selection.chr2L.vs_droSim1.bwaUniq.vcf --out selection.chr2L.vs_droSim1.bwaUniq --freq 


#bedtools intersect -wa -wb -a <(cat PopSech.chr2L.vs_droSim1.bwaUniq.frq | tail -n +2 | awk '{print $1,$2,$2+1,$4,$5,$6}' | tr " " "\t") -b <(cat PopSim.chr2L.vs_droSim1.bwaUniq.frq | tail -n +2 | awk '{print $1,$2,$2+1,$4,$5,$6}' | tr " " "\t") | bedtools intersect -wa -wb -a - -b <(cat selection.chr2L.vs_droSim1.bwaUniq.frq | tail -n +2 | awk '{print $1,$2,$2+1,$4,$5,$6}' | tr " " "\t") | cut -f 1,2,4-6,10,12,16,18 | tr ":" "\t" | awk '{print $1,$2,$2+1,"0","0","+",$4,$6,$3,$7,$8,$10,$11,$13}' | tr " " "\t" > dev/grouped_freaks2.bed


#*.bwaUniq.frq ---> grouped_freaks2

#( dev/grouped_freaks2.bed --> R script --> chr2L.freqCompare.bed)

#cat chr2L.freqCompare.bed | cut -f 1,2,3,6,18,19 | nl | tr -d " " | awk '{print $2,$3,$4,$1,"0",".",$6,$7}' | tr " " "\t" > chr2L.freqCompare.clean.bed 

#bedtools map -c 7,8,8 -o sum,sum,count -null NA -a dev/droSim1_w100000_s100000.windows.chr2L.bed -b chr2L.freqCompare.clean.bed > dev/chr2L.freqCompare.windowed.bed

############################################################################  




#######		Calculating VCF-based distance metric by windowed 	####

rule pairwise_VCF_distance_metric_windowed:
	input:
		vcf_in="variants/{prefix}.vs_{ref_genome}.{aligner}.vcf",
		windoze_in = "utils/{window_prefix}.windows.bed" #generalize this
	output:
		pairdist_out = ["variants/{prefix}.vs_{ref_genome}.{aligner}/distances/{indiv_1}/{indiv_2}.{window_prefix}.bed","variants/{prefix}.vs_{ref_genome}.{aligner}/distances/{indiv_2}/{indiv_1}.{window_prefix}.bed"]
	params:
		runmem_gb=16,
		runtime = "1:00:00",
	run:
		shell(
			"""
			mkdir -p variants/{wildcards.prefix}.vs_{wildcards.ref_genome}.{wildcards.aligner}/distances/{wildcards.indiv_1}/
			bedtools map -a {input.windoze_in} -b <(vcf-subset {input.vcf_in} -u -c {wildcards.indiv_1},{wildcards.indiv_2} | grep -v "#" | sed -e 's/0\/0:\S*\(\s\|$\)/0\t/g'  | sed -e 's/0\/1:\S*\(\s\|$\)/1\t/g' | sed -e 's/1\/0:\S*\(\s\|$\)/1\t/g' | sed -e 's/1\/1:\S*\(\s\|$\)/2\t/g' | sed -e 's/\.:\S*\(\s\|$\)/NA\t/g' |awk '{{print $1,$2,$2+1,0,0,"*",$10,$11,sqrt(($10-$11)^2)}}' | tr " " "\t" | grep -v NA ) -c 9,9 -o sum,count | awk '{{if($6>0)print $1,$2,$3,$4,$5/(2*$6) ;else print $1,$2,$3,$4,"NA"}}' | tr " " "\t" > {output.pairdist_out[0]}
			"""
			)
		if len(set(output.pairdist_out)) > 1:
			shell(
				"""
				mkdir -p variants/{wildcards.prefix}.vs_{wildcards.ref_genome}.{wildcards.aligner}/distances/{wildcards.indiv_2}/
				cp {output.pairdist_out[0]} {output.pairdist_out[1]}
				"""
				)
# 			hmm..... maybe wrap that one-liner up in a different script?
#all_samples.chr2L.vs_droSim1.bwaUniq.vcf

rule all_dist_to_indiv_by_group:
	input:
		groupies = lambda wildcards: expand("variants/{prefix}.vs_{ref_genome}.{aligner}/distances/{indiv}/{groupie}.{window_prefix}.bed", prefix=wildcards.prefix, ref_genome=wildcards.ref_genome, aligner=wildcards.aligner, indiv=wildcards.indiv, groupie=sampname_by_group[wildcards.group], window_prefix=wildcards.window_prefix),
	output:
		group_dist = "variants/{prefix}.vs_{ref_genome}.{aligner}/distances/{indiv}/{group}.{window_prefix}.tbl"
	params:
		runmem_gb=8,
		runtime = "10:00",
	run:
		shell("touch {output.group_dist}.tmp")

		for groupie in sampname_by_group[wildcards.group]:
			shell("""cat variants/{wildcards.prefix}.vs_{wildcards.ref_genome}.{wildcards.aligner}/distances/{wildcards.indiv}/{groupie}.{wildcards.window_prefix}.bed | cut -f 5 | cat <(echo {groupie}) - | paste {output.group_dist}.tmp - > {output.group_dist} 
					cp {output.group_dist} {output.group_dist}.tmp
					""")
		shell("""cat <(echo "chr\tstart\tstop\twin") variants/{wildcards.prefix}.vs_{wildcards.ref_genome}.{wildcards.aligner}/distances/{wildcards.indiv}/{groupie}.{wildcards.window_prefix}.bed | cut -f 1-4 | paste - {output.group_dist}.tmp > {output.group_dist}""")
		shell("rm {output.group_dist}.tmp")

rule group2group_VCF_distance:
	input:
		tables_in = lambda wildcards: expand("variants/{prefix}.vs_{ref_genome}.{aligner}/distances/{indiv}/{group2}.{window_prefix}.tbl", prefix=wildcards.prefix, ref_genome=wildcards.ref_genome, aligner=wildcards.aligner, group2=wildcards.group2, indiv=sampname_by_group[wildcards.group1], window_prefix=wildcards.window_prefix)
	output:
		roster_out = "variants/{prefix}.vs_{ref_genome}.{aligner}/distances/from.{group1}.to.{group2}.{window_prefix}.distanceRoster"
	params:
		runmem_gb=1,
		runtime = "1:00",
	shell:
		"""echo {input.tables_in} | tr " " "\n" > {output.roster_out}"""

############################################################


rule write_report:
	input:
		reference_genome_summary = ["meta/reference_genomes.summary"],
		sequenced_reads_summary=["meta/sequenced_reads.dat"],
		alignment_summaries = expand("meta/alignments.vs_{ref_genome}.{aligner}.summary", ref_genome=['droSim1', 'droSec1'], aligner=['bwa','bwaUniq']),
		full_variant_summary = expand("meta/{prefix}.calledVariants.{aligner}.summary", aligner=["bwaUniq"], prefix=["all_samples"] ),
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






