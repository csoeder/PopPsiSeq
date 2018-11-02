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

rule dummy:
	input:
		filesin = return_filename_by_sampname,
	output:
		header="{samplename}.head"
	shell:
		"head {input.filesin} > {output.header}"
	


