configfile: 'config.yaml'

sample_by_name = {c['name'] : c for c in config['data_sets']}


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
		reportBack = "meta/FASTP/{samplename}.json"
	params:
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json meta/FASTP/{samplename}.json --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in1 {input.fileIn[1]} --out1 {output.fileOut[1]}"

#request reportBacks instead, process them into a unified SQL upload?
rule clean_all_samps:
	input: 
		clean_se = [sample_by_name[nom]['path']+nom+".clean.R0.fastq" for nom in sample_by_name.keys() if not sample_by_name[nom]['paired']],
		clean_pe = [sample_by_name[nom]['path']+nom+".clean.R"+arr+".fastq" for nom in sample_by_name.keys() if sample_by_name[nom]['paired'] for arr in ["1","2"]],
	output:
		outflg = "allclean.flag"
	shell:
		"touch {output.outflg}"
