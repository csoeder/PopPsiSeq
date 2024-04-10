import argparse
import re
from numpy import argmax

parser = argparse.ArgumentParser()
parser.add_argument("hybrid_pileup", help="pileup file for hybrid")
parser.add_argument("parent_pileup", help="pileup file for parent")
parser.add_argument("output", help="file to write shared SNPs to")
parser.add_argument("-v", "--verbose", action="count", default=0, help="verbose reporting")
parser.add_argument("-c", "--minCov", help="minimum coverage for a variant to be considered", default=3, type=int)
parser.add_argument("-f", "--minFrac", help="minimum agreement fraction for a variant to be considered", default=0.9, type=float)
parser.add_argument("-F", "--fai_file", help=".fai file defining sort order")

args = parser.parse_args()

verbose = args.verbose
minCov = args.minCov
minFrac = args.minFrac
fai_file = args.fai_file
indel_regex = re.compile('[\+\-][0-9]+[ACGTNacgtn]+')
mismatch_regex = re.compile('[ACGTacgt]')


def get_sort_order(fai):
	sorta = []
	phial = open(fai, "r")
	for line in (phial):
		sorta.append(line.split("\t")[0])
	return sorta


def scan_next_contig(file_in):
	
	global sort_order

	# if this is the first line, start there
	# otherwise start where we left off
	if file_in.name not in bookmarks.keys():
		lion = file_in.readline()
	else:
		lion = bookmarks[file_in.name]

	next_contig = lion.split("\t")[0]

	if next_contig != "":
		#if we haven't reached the EOF

		#if this contig was eliminated from the sort_order by emptiness in the other pileup, skip it
		while next_contig not in sort_order:
			lion = file_in.readline()
			next_contig = lion.split("\t")[0]

		#is this the next in the sort order? probably, but if not, then any before this point are irrelevent. skip them
		# ndx = sort_order.index(next_contig)
		# if verbose > 2 and ndx > 0:
			# print("skipping %s empty contigs" % (ndx) )
		# sort_order = sort_order[ndx:]

		if verbose > 1:
			print("scanning contig %s in pileup %s" % (next_contig, file_in.name) )

		contig_dict = dict((
											('mini_counter', 0),
											('sum_coverage', 0),
											('mm_count', 0),
											('mini_counter', 0),
											('position_dict', {}),
											('contig', next_contig)))

		while lion.split("\t")[0] == next_contig:
			


			try:
				contig, position, ref_base, coverage, read_base, qual = lion.split()
			except ValueError:  # zero-cov sites don't report a readbase/qual Dx
				contig, position, ref_base, coverage = lion.split()
				read_base, qual = "N", "I"


			contig_dict['mini_counter'] += 1
			contig_dict['sum_coverage'] += int(coverage)


			if (
				int(coverage) > 1 and "*" not in read_base and not
				indel_regex.search(read_base) and ref_base.upper() in ['A','T','C','G'] and
				mismatch_regex.search(read_base)):
				# if cov >1 and read_base is not an indel and the reference base
				# is not null and read bases contain mismatches
				# ... then tally the bases at this site
				read_base_dict = {
									"A": read_base.upper().count("A"),
									"T": read_base.upper().count("T"),
									"C": read_base.upper().count("C"),
									"G": read_base.upper().count("G")}
				read_base_dict[ref_base.upper()] += (
													read_base.upper().count(",") +
													read_base.upper().count("."))
				#  incremenet when the read base is the reference base
				contig_dict['mm_count'] += 1  # increment mismatch count
				contig_dict['position_dict'][int(position)] = {'ref_base':ref_base.upper(), 'pileup_cov':int(coverage), 'base_dict':read_base_dict} 

			if int(position) % 100000 == 0 and verbose > 0:
				print("Contig %s: %s KiloBases scanned!" % tuple([contig, int(position)/1000]))

			lion = file_in.readline()

		bookmarks[file_in.name] = lion
		return contig_dict
	else: # if we're at the EOF
		return None


def contig_report(contig_dict):
	if contig_dict['mm_count'] > 0:
		mismatch_warn_string = '.'
	else:
		mismatch_warn_string = ' (only Ns in the pileup reference base column?)'
	print( "contig %s had an average coverage depth of %s reads and a raw mismatch count of %s%s" % tuple([contig_dict["contig"], contig_dict['sum_coverage']/contig_dict['mini_counter'],  contig_dict['mm_count'], mismatch_warn_string]) )

	
def mismatch_chooser(site_dict):
	#{'base_dict': {
	#	'A': 0, 'C': 0, 'T': 0, 'G': 28}, 'ref_base': 'A', 'pileup_cov': 28}}
	# choice = 'N'
	choice = site_dict['ref_base']
	meta = None # 	We may want to later report information on polymorphism
	for bass in site_dict['base_dict'].keys():
		if int(float(site_dict['pileup_cov'])) >= minCov and site_dict['base_dict'][bass]/float(site_dict['pileup_cov']) >= minFrac :
			choice = bass
	return choice, meta


def contig_comparator(parent_dict, hybrid_dict):

	comparisons = []

	if verbose > 0:
		print("starting the parent-hybrid comparison comparison for contig %s" % (parent_dict["contig"]) )

	minicount = 0
	total_count = len(parent_dict['position_dict'].keys())

	if verbose > 0:
		print("finding common positions...")

	common = set(parent_dict['position_dict'].keys()).intersection(set(hybrid_dict['position_dict'].keys()))
	if verbose > 0:

		print("common positions found!!")


	for parent_pos in common:

		minicount += 1			

		parent_minidict = parent_dict['position_dict'][parent_pos]

		par_var, par_meta = mismatch_chooser(parent_minidict)

		if par_var != parent_minidict['ref_base']:
				#if the site isn't actually variable in the parent, ignore it
			hybrid_minidict = hybrid_dict['position_dict'][parent_pos]

			hyb_var, hyb_meta = mismatch_chooser(hybrid_minidict)

			if hybrid_minidict['ref_base'] !=  parent_minidict['ref_base']:
				# If this happens... something, somewhere has gone terribly wrong x_x
				print( "WARNING: reference sequences disagree on contig %s, position %s !!!" % tuple([contig, parent_pos]))

			elif hyb_var == par_var:
					#if the site has the same variant in hybrid and parent, record it as parent-derived
				comparisons.append([parent_pos, 1])
			else:
				#if the parent site is the wrong variant in the hybrid, it's not parent-derived
				comparisons.append([parent_pos, 0])
#		else: #?????????????????????????????????????????????????????????
		# 	If the parent variant site isn't variant in the hybrid, the hybrid site isn't parent-derived.
#			comparisons.append([parent_pos, 0])
			#parent_minidict = parent_dict['position_dict'].pop(parent_pos)


		if minicount % 1000 == 0 and verbose > 0:
			print( "%s parent mismatch sites investigated of %s!" % tuple([minicount, total_count]))

	return comparisons


def comparison_writer(comps, contig, file_out):
	for coord in comps:
		file_out.write("%s\t%s\t%s\t%s\n" % tuple([contig, coord[0], coord[0]+1, coord[1]]) )






bookmarks = {}
sort_order = get_sort_order(fai_file)
parent_mpile = open(args.parent_pileup, "r")
hybrid_mpile = open(args.hybrid_pileup, "r")
write_file = open(args.output, "w")

piles = [ parent_mpile, hybrid_mpile ]

dix = [ scan_next_contig(parent_mpile), scan_next_contig(hybrid_mpile) ]

while None not in dix:

	while dix[0]["contig"] != dix[1]["contig"] :
		advanced = argmax([sort_order.index(dix[0]["contig"]), sort_order.index(dix[1]["contig"])  ] ) 
		# which one is farther along the sort order, ie has skipped the most?

		dix[1-advanced] = scan_next_contig(piles[1-advanced])
		# move the other one up 

	if verbose >0:
		print("done scanning contig %s in %s!" % (dix[0]["contig"], piles[0].name))
		contig_report(dix[0])
		print()
		print("done scanning contig %s in %s!" % (dix[1]["contig"], piles[1].name))
		contig_report(dix[1])
		print()
		print()
		print()


	comparisons = contig_comparator(dix[0], dix[1])

	comparison_writer( comparisons , dix[0]["contig"] ,  write_file)

	dix = [ scan_next_contig(parent_mpile), scan_next_contig(hybrid_mpile) ]


parent_mpile.close()
hybrid_mpile.close()
write_file.close()


print("and that's all the contigs!")
	



# sort_order = get_sort_order("/proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_sechellia/droSec1.fa.fai")
# parent_mpile = open("test1.mpile", "r")
# hybrid_mpile = open("test2.mpile", "r")

# bookmarks = {}
# sync = 0 
# piles = [ parent_mpile, hybrid_mpile ]






print( "DONE!")
