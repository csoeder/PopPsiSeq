#configfile: 'config.yaml'

from Bio import SeqIO

# module load python/3.9.6 sratoolkit samtools freebayes vcftools bwa bedtools r/4.2.2 rstudio
#PATH=$PATH:/nas/longleaf/home/csoeder/modules/vcflib/bin:/nas/longleaf/home/csoeder/modules/parallel/bin


sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}

try:
	chain_dict_by_destination = config['lift_genomes']
except KeyError:
	chain_dict_by_destination = {}

sampname_by_group = {}
sampname_by_sra = {}
#print(sample_by_name)

for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]

	try:
		sampname_by_sra[sample_by_name[s]['SRA']] = s
	except KeyError:
		pass


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

############################################################################  
#######		background data, eg reference genomes 	####
############################################################################  


rule reference_genome_reporter:
	input:
		fai_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['fai'],
	output:
		report_out = "data/summaries/reference_genomes/{ref_gen}.fai.report"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"""
		mkdir -p data/summaries/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule demand_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("data/summaries/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "data/summaries/reference_genomes/reference_genomes.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"cat {input.refgen_reports} > {output.refgen_summary}"

############################################################################  
#######		Read files: summon them, process them 	####
############################################################################  

rule summon_reads_SRA_pe:
	output:
		reads1='data/external/sequence/paired_end/{prefix}/{prefix}_1.fastq',
		reads2='data/external/sequence/paired_end/{prefix}/{prefix}_2.fastq',
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
	run:
		try:
			sra = sampname_by_sra[wildcards.prefix]#["SRA"]
			shell(""" mkdir -p data/external/sequence/paired_end/{wildcards.prefix}/ """)
			shell("""
				fasterq-dump  --split-3 --outdir data/external/sequence/paired_end/{wildcards.prefix}/ {wildcards.prefix}
			""")
		except KeyError:
			raise KeyError("Sorry buddy, you can only download SRAs that are associated with a sample in the config file! " )

rule summon_reads_SRA_se:
	output:
		reads='data/external/sequence/single_end/{prefix}/{prefix}.fastq',
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
	run:

		try:
			sra = sampname_by_sra[wildcards.prefix]#["SRA"]
			shell(""" mkdir -p data/external/sequence/single_end/{wildcards.prefix}/ """)
			# shell("""
				# fasterq-dump  --split-3 --outdir data/external/sequence/single_end/{wildcards.prefix}/ {sra}
			# """)
			shell("""
				fasterq-dump  --split-3 --outdir data/external/sequence/single_end/{wildcards.prefix}/ {wildcards.prefix}
			""")

		except KeyError:
			raise KeyError("Sorry buddy, you can only download SRAs that are associated with a sample in the config file! " )


def return_file_relpath_by_sampname(sampname):


	if sample_by_name[sampname]["source"] in ["NCBI"]:
		subfolder = sample_by_name[sampname]["SRA"]
		path = "data/external/sequence/"

	else:
		subfolder = sampname
		path = "data/raw/sequence/"

	if sample_by_name[sampname]["paired"] :
		path = "%spaired_end/%s/" % (path, subfolder)

	else:
		path = "%ssingle_end/%s/" % (path,subfolder)

	filesin = return_filename_by_sampname(sampname)
	pathsout = ["".join([path, fq]) for fq in filesin]

	return pathsout

def return_file_relpath_by_sampname(sampname):


	pith = sample_by_name[sampname]["path"]

	filesin = return_filename_by_sampname(sampname)
	pathsout = ["".join([pith, fq]) for fq in filesin]

	return pathsout


rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards.samplename)
	output:
		fileOut = ["data/intermediate/sequence/{samplename}/{samplename}.clean.R0.fastq"],
#		fileOut = ["{pathprefix}/{samplename}.clean.R0.fastq"],
		jason = "data/intermediate/sequence/{samplename}/{samplename}.False.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json data/intermediate/sequence/{samplename}/{samplename}.False.json --n_base_limit 0 ",# --html meta/FASTP/{samplename}.html", 
		se_params = "",
	message:
		"FASTP QA/QC on single-ended reads ({wildcards.samplename}) in progress.... "
	run:
		shell(""" mkdir -p data/intermediate/sequence/{wildcards.samplename}/ """)
		shell(""" /nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} """)
		
rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards.samplename)
	output:
		fileOut = ["data/intermediate/sequence/{samplename}/{samplename}.clean.R1.fastq","data/intermediate/sequence/{samplename}/{samplename}.clean.R2.fastq"],
		jason = "data/intermediate/sequence/{samplename}/{samplename}.True.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json data/intermediate/sequence/{samplename}/{samplename}.True.json --n_base_limit 0 ",# --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	message:
		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
	run:
		shell(""" mkdir -p data/intermediate/sequence/{wildcards.samplename}/ """)
		shell(""" /nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]} """)


rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("data/intermediate/sequence/{samp}/{samp}.{pairt}.json", samp = wildcards.samplename, pairt = sample_by_name[wildcards.samplename]['paired'])
#		jason = lambda wildcards: expand("{path}{samp}.{pairt}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, pairt = sample_by_name[wildcards.samplename]['paired'])
	output:
		jason_pruned = "data/summaries/intermediate/FASTP/{samplename}/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		mkdir -p data/summaries/intermediate/FASTP/{wildcards.samplename}/
		cp {input.jason} data/summaries/intermediate/FASTP/{wildcards.samplename}/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""

rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = lambda wildcards: expand("data/summaries/intermediate/FASTP/{samplename}/{samplename}.json.pruned", samplename = sampname_by_group[wildcards.group])
	output:
		summary = "data/summaries/intermediate/FASTP/{group}.sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"



############################################################################  
#######		Map reads to reference genomes 	####
############################################################################  


rule bwa_align:
	input:
#		reads_in = lambda wildcards: expand("data/intermediate/FASTQs/{source}/{sample}/{sample}.clean.R{arr}.fastq", source=sample_by_name[wildcards.sample]['source'], sample=wildcards.sample, arr=[ [1,2] if sample_by_name[wildcards.sample]['paired'] else [0] ][0]),
		reads_in = lambda wildcards: expand("data/intermediate/sequence/{sample}/{sample}.clean.R{arr}.fastq", sample=wildcards.sample, arr=[ [1,2] if sample_by_name[wildcards.sample]['paired'] else [0] ][0]),
		ref_genome_file = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]['path'],
	output:
		bam_out = "data/intermediate/mapped_reads/bwa/{sample}.vs_{ref_genome}.bwa.sort.bam",
	params:
		runmem_gb=96,
		runtime="64:00:00",
		cores=8,
	message:
		"aligning reads from {wildcards.sample} to reference_genome {wildcards.ref_genome} .... "
	run:
		shell("bwa aln {input.ref_genome_file} {input.reads_in[0]} > {input.reads_in[0]}.{wildcards.ref_genome}.sai ")
		if sample_by_name[wildcards.sample]['paired']:
			shell("bwa aln {input.ref_genome_file} {input.reads_in[1]} > {input.reads_in[1]}.{wildcards.ref_genome}.sai ")
			shell("bwa sampe {input.ref_genome_file} {input.reads_in[0]}.{wildcards.ref_genome}.sai {input.reads_in[1]}.{wildcards.ref_genome}.sai {input.reads_in[0]}  {input.reads_in[1]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		else:
			shell("bwa samse {input.ref_genome_file} {input.reads_in[0]}.{wildcards.ref_genome}.sai {input.reads_in[0]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		shell("samtools index {output.bam_out}")

rule bwa_uniq:
	input:
		bam_in = "data/intermediate/mapped_reads/bwa/{sample}.vs_{ref_genome}.bwa.sort.bam"
	output:
		bam_out = "data/intermediate/mapped_reads/bwaUniq/{sample}.vs_{ref_genome}.bwaUniq.sort.bam"
	params:
		quality="-q 20 -F 0x0100 -F 0x0200 -F 0x0300 -F 0x04",
		uniqueness="XT:A:U.*X0:i:1.*X1:i:0",
		runmem_gb=16,
		runtime="18:00:00",
		cores=4,
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
		bam_in = "data/intermediate/mapped_reads/{aligner}/{sample}.vs_{ref_genome}.{aligner}.sort.bam"
	output:
		report_out = "data/summaries/intermediate/BAMs/{aligner}/{sample}.vs_{ref_genome}.{aligner}.summary",
		dpth_by_chrom = "data/intermediate/mapped_reads/{aligner}/{sample}.vs_{ref_genome}.{aligner}.sort.bam.dpth_by_chrom",
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=1,
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	run:
		shell(""" mkdir -p data/summaries/intermediate/BAMs/{wildcards.aligner}/ """)
		ref_genome_idx=ref_genome_by_name[wildcards.ref_genome]['fai']
		shell("samtools idxstats {input.bam_in} > {input.bam_in}.idxstats")
		shell("samtools flagstat {input.bam_in} > {input.bam_in}.flagstat")
		shell("bedtools genomecov -max 1 -ibam {input.bam_in} -g {ref_genome_idx} > {input.bam_in}.genomcov")
		#change the -max flag as needed to set 
		shell("""samtools depth -a {input.bam_in} | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' > {input.bam_in}.dpthStats""")


		shell("""samtools depth -a {input.bam_in}  > {input.bam_in}.dpth.tmp """)
                shell(""" 
                rm -rf {input.bam_in}.dpth_by_chrom;
                for chrom in $(cat {ref_genome_idx} | cut -f 1); do 
                cat {input.bam_in}.dpth.tmp | grep -w "$chrom" | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk -v chr="$chrom" '{{print"{wildcards.sample}\t"chr"\t"$0}}' >>  {input.bam_in}.dpth_by_chrom ; 
                done  """)

                shell("""rm {input.bam_in}.dpth.tmp """)



		#https://www.biostars.org/p/5165/
		#save the depth file and offload the statistics to the bam_summarizer script?? 
		shell("python3 scripts/bam_summarizer.py -f {input.bam_in}.flagstat -i {input.bam_in}.idxstats -g {input.bam_in}.genomcov -d {input.bam_in}.dpthStats -o {output.report_out} -t {wildcards.sample}")




rule demand_BAM_analytics:
	input:
		bam_reports = lambda wildcards: expand("data/summaries/intermediate/BAMs/{aligner}/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sampname_by_group[wildcards.group], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
                dpth_reports = lambda wildcards: expand("data/intermediate/mapped_reads/{aligner}/{sample}.vs_{ref_genome}.{aligner}.sort.bam.dpth_by_chrom", sample=sampname_by_group[wildcards.group], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),

	output:
#		full_report = "data/summaries/intermediate/BAMs/{group}.vs_{ref_genome}.{aligner}.summary"
		full_report = "data/summaries/intermediate/BAMs/{group, \w+}.vs_{ref_genome}.{aligner}.summary",
		all_dpth_by_chrom = "data/summaries/intermediate/BAMs/{group, \w+}.vs_{ref_genome}.{aligner}.dpth_by_chrom"

	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"cat {input.bam_reports} > {output.full_report}; cat {input.dpth_reports} > {output.all_dpth_by_chrom}"


############################################################################  
#######		Call variants from mapped reads 	####
############################################################################  


rule joint_vcf_caller_parallel:
	input:
		bams_in = lambda wildcards: expand("data/intermediate/mapped_reads/{aligner}/{sample}.vs_{ref_genome}.{aligner}.sort.bam", sample=sampname_by_group[wildcards.group], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
		windows_in = "utils/genome_windows/{ref_genome}_w100000_s100000.windows.list"
	output:
		vcf_out = "data/intermediate/variants/freebayes/{group}.vs_{ref_genome}.{aligner}.vcf"
	params:
		freebayes="--standard-filters",
		runmem_gb=64,
		runtime="10-00:00:00",
		cores=24,
	message:
		"Using Freebayes to jointly call variants from all samples mapped to \ {wildcards.ref_genome} \ with \ {wildcards.aligner} \ "
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		
		shell("scripts/freebayes-parallel {input.windows_in} {params.cores} {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out} ")


rule vcf_reporter: # replace this with the better reporter for freebayesSmrt?
	input:
		vcf_in = "data/intermediate/variants/{prefix}.vs_{ref_genome}.{aligner}.vcf"
	output:
		report_out = "data/summaries/intermediate/VCFs/{ref_genome}/{prefix}.vs_{ref_genome}.{aligner}.summary",
		frq_out = "data/intermediate/frequencies/{prefix}.vs_{ref_genome}.{aligner}.frq"
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=4,
#	message:
#		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	shell:
		"""
		mkdir -p data/intermediate/frequencies/ data/summaries/intermediate/VCFs/{wildcards.ref_genome}/
		cat {input.vcf_in}  | grep -v "#" | cut -f 1 | sort | uniq -c > {output.report_out}.snpsPerContig.tmp
		cat {output.report_out}.snpsPerContig.tmp | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' | cat {output.report_out}.snpsPerContig.tmp - > {output.report_out}.snpsPerContig
		rm {output.report_out}.snpsPerContig.tmp
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --freq 
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --counts
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-indv
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-site
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --singletons

		cat {output.report_out}.frq | tail -n +2 | awk '{{print $1,$2,$2+1,$4,$5,$6}}' | tr " " "\t" > {output.frq_out}

		ref_genome={wildcards.ref_genome}

		tail -n 1 {output.report_out}.snpsPerContig | awk '{{print "total_snp_count\t"$1}}' | sed -e 's/^/'$ref_genome'\t/g' > {output.report_out}
		"""


rule summon_VCF_analytics_base:
	input:
		bam_reports = lambda wildcards: expand("data/summaries/intermediate/VCFs/{ref_genome}/{caller}/{prefix}.vs_{ref_genome}.{aligner}.summary", prefix=wildcards.prefix, ref_genome=["dm6","droSim1","irvineSec"], aligner="bwaUniq")
	output:
		full_report = "data/summaries/intermediate/VCFs/{caller}/{prefix}.calledVariants.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"""
		prefix={wildcards.prefix}
		cat {input.bam_reports} | sed -e 's/^/'$prefix'\t/g'> {output.full_report}
		"""


rule subset_VCF_to_subgroup: # does this work with freebayesSmrt?
	input:
		vcf_in = "data/intermediate/variants/{caller}/{group}.vs_{ref_genome}.{aligner}.vcf"
	output:
		vcf_out = "data/intermediate/variants/{caller}/{group}.subset_{subgroup}.vs_{ref_genome}.{aligner}.vcf"
	params:
		runmem_gb=8,
		runtime = "1:00:00",
		cores=2,
	message:
		"Subsetting the variant file \ {input.vcf_in} \ to the individuals in the subgroup \ {wildcards.subgroup} \ .... "
	run:
#		filter_string = "--min-alleles 2 --max-alleles 2  --max-missing-count 0 " #--non-ref-ac-any 1" # only biallelic sites, no missing genotypes, no invariant sites.
#		member_list = "%s,"*len(sampname_by_group[wildcards.subgroup]) % tuple(sampname_by_group[wildcards.subgroup])
#		shell("vcf-subset {input.vcf_in} -u -c {member_list} | vcftools {filter_string} --vcf - --recode --stdout > {output.vcf_out}")

		member_list = "%s,"*len(sampname_by_group[wildcards.subgroup]) % tuple(sampname_by_group[wildcards.subgroup])
		filter_string = " --min-alleles 2 --max-alleles 2 --types snps --exclude 'F_MISSING>0' "
		shell("vcf-subset {input.vcf_in} -u -c {member_list} | bcftools +fill-tags -  -- -t AF,AC,F_MISSING,NS,AN | bcftools view {filter_string}  > {output.vcf_out}")






rule window_maker:
	output:
		windowed='utils/genome_windows/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed',
		winlist='utils/genome_windows/{ref_genome}_w{window_size}_s{slide_rate}.windows.list',

	params:
		runmem_gb=8,
		runtime="5:00",
		cores=1,
	run:
		fai_path = ref_genome_by_name[wildcards.ref_genome]['fai'],
		shell("mkdir -p utils/genome_windows/")
		shell(
			'bedtools makewindows -w {wildcards.window_size} -s {wildcards.slide_rate} -g {fai_path} -i winnum | bedtools sort -i - > {output.windowed}'
		)
		shell("""cat {output.windowed}| awk '{{print$1":"$2"-"$3}}' > {output.winlist}""")

############################################################################  
#######		windowed frequency comparison between groups 	################
############################################################################  

rule calc_frq_shift:
	input:
		par1_frq = "data/intermediate/frequencies/{prefix}.subset_{grup_par1}.vs_{ref_genome}.{aligner}.frq",
		par2_frq = "data/intermediate/frequencies/{prefix}.subset_{grup_par2}.vs_{ref_genome}.{aligner}.frq",
		off_frq = "data/intermediate/frequencies/{prefix}.subset_{grup_off}.vs_{ref_genome}.{aligner}.frq",
	output:
		frqShft_out = "data/intermediate/freq_shift/{prefix}.{grup_off}_with_{grup_par1}_and_{grup_par2}.vs_{ref_genome}.{aligner}.frqShift"
	params:
		runmem_gb=16,
		runtime="1:00:00",
		cores=2,
	shell:
		"""
		mkdir -p data/intermediate/freq_shift/
		bedtools intersect -wa -wb -a {input.par1_frq} -b  {input.par2_frq} | bedtools intersect -wa -wb -a - -b  {input.off_frq} | cut -f 1,2,4-6,10,12,16,18 | tr ":" "\t" | awk '{{print $1,$2,$2+1,"0","0","+",$4,$6,$3,$7,$8,$10,$11,$13}}' | tr " " "\t" > {output.frqShft_out}.pre
		Rscript scripts/freqShifter.R {output.frqShft_out}.pre {output.frqShft_out}
		rm {output.frqShft_out}.pre
		"""



rule window_frq_shift: #generalized windows? eg from moehring paper
	input:
		frqShft_in = "data/intermediate/freq_shift/{frqshft_prefix}.vs_{ref_genome}.{aligner}.frqShift",
		windows_in = "utils/genome_windows/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed"
	output:
		windowed_out = "data/ultimate/freq_shift/{frqshft_prefix}.vs_{ref_genome,[^.]+}.{aligner}.windowed_w{window_size}_s{slide_rate}.frqShift.bed"
	params:
		runmem_gb=16,
		runtime="10:00",
		cores=4,
	shell:
		"""
		mkdir -p data/ultimate/freq_shift/
		bedtools map -c 7,8,8 -o sum,sum,count -null NA -a {input.windows_in} -b <( tail -n +2  {input.frqShft_in} | cut -f  1-3,15,16 | nl | tr -d " " | awk '{{if( $5!="NA" && $6!="NA")print $2,$3,$4,$1,"0",".",$5,$6}}' | tr " " "\t"| bedtools sort -i - ) > {output.windowed_out}
  		"""

rule bin_frqShft:
	input:
		shared_snps='data/intermediate/freq_shift/{var_caller}/{frqshft_prefix}.frqShift',
	output:
		windowed_out = "data/ultimate/freq_shift/{var_caller}/{frqshft_prefix}.snpBinned_b{bin_size}_s{slide_rate}.frqShift.bed"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	run:
		shell("python3 scripts/bin_by_SNPs.py -i {input.shared_snps} -o {output.windowed_out}.tmp -b {wildcards.bin_size} -s {wildcards.slide_rate}  -c 15,16,16 -m 'sum,sum,count' --header_skip 1 " )
		shell("paste <( cut -f 1-3 {output.windowed_out}.tmp ) <( cut -f 4- {output.windowed_out}.tmp | nl -n rz ) > {output.windowed_out}")
		shell("rm {output.windowed_out}.tmp " )



rule PopPsiSeq_lifter:
	input:
#		unlifted='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.sharedSnps.bed'
		unlifted="data/intermediate/freq_shift/{prefix}.vs_{ref_genome}.{aligner}.frqShift"
	output:
#		lifted='data/intermediate/freq_shift/PsiSeq{version,.*}/{lift_genome}/{aligner}/{sample}.SNPs_shared_with.{eff_naught}.vs_{ref_genome}.{aligner}.lifted_to_{lift_genome}.bed',
		lifted='data/intermediate/freq_shift/{prefix}.vs_{lift_genome}.{aligner}.lifted_from_{ref_genome}.frqShift',
		too_heavy='data/intermediate/freq_shift/{prefix}.vs_{lift_genome}.{aligner}.2Heavy2Lift_from_{ref_genome}'
	params:
		runmem_gb=32,
		runtime="12:00:00",
		cores=1
	run:
#		chain = chain_dict_by_destination[wildcards.ref_genome][wildcards.lift_genome]
		chain = chain_dict_by_destination[wildcards.lift_genome][wildcards.ref_genome]
#		shell("mkdir -p data/intermediate/shared_snps/PsiSeq{wildcards.version}/{wildcards.lift_genome}/{wildcards.aligner}/")
		shell(
			'/nas/longleaf/home/csoeder/modules/UCSC_utils/liftOver -bedPlus=3 <( tail -n +2 {input.unlifted} ) {chain} {output.lifted}.tmp {output.too_heavy}'
		)
		shell(
			'bedtools sort -i {output.lifted}.tmp > {output.lifted}'
		)
		shell(
			'rm {output.lifted}.tmp'
		)

ruleorder: calc_frq_shift > PopPsiSeq_lifter





rule DAGnabbit:
	input:
		target = "{path_prefix}/{file_prefix}"
	output:
		dag_nabbed = "meta/pipeline_schematics/{path_prefix}/{file_prefix,[^\/]+}.png"
	params:
		runmem_gb=1,
		runtime="10:00",
		cores=8,
	run:
		shell("""rm -rf {output.dag_nabbed};  mkdir -p meta/pipeline_schematics/{wildcards.path_prefix}/ """)

		shell(""" snakemake --rulegraph {input.target} | dot -T png > {output.dag_nabbed} """)


