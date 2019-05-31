"""
Owner: Ada Madejska, UCSB, J.Rothman's lab 2019
Description:
	This program can be used for looking for SNPs frequency across multiple organisms
	of the same species. It creates a graph of specified region from the genome and 
	looks at the frequency at which the nucleotides differ from the reference at each
	position. 
Input:
	Sorted and indexed BAM file (1) with all sequenced and aligned genomes for study
	Sorted and indexed FASTA reference file with which the BAM files were aligned.
	GFF3 gzipped gene annotation file if the user wants to find gene names for mutations. 
Output:
	A graph showing the fractions of mutations happening at each position in the 
	specified regions.
	List of genes affected by the most frequent SNPs with number of mutations.
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import sys


def find_saturation(bam, ref, start, end, chrom, rs, re):
	"""
	Reads the BAM file and counts each base at a specific aligned position.
	Compares those reads to the reference and calculates frequency of the SNPs at each
	position.
	:param string bam: BAM file pathway.
	:param string ref: Reference file pathway.
	:param int start: Start position.
	:param int end: End position. 
	:param float rs: Range start number.
	:param float re: Range end number.
	:return: a dictionary of mutations at each position.
	"""

	# Read the BAM file.
	bamfile = pysam.AlignmentFile(bam, 'rb')

	# Read reference FASTA file.
	fastafile = pysam.FastaFile(ref)

	mutations = {}
	total_reads = 0

	# fetch() returns all reads overlapping a region sorted by the first aligned base in the reference
	# 	sequence. Note that it will also return reads that are only partially overlapping with 
	#	the region.
	# Create the dictionary of dictionaries for SNPs at each position.
	for read in bamfile.fetch(chrom, start, end):

		# read.positions gives an array of all the positions of each sequence.
		positions = read.get_reference_positions()
		sequence = read.query_alignment_sequence # Don't want soft clipped bases.
		quality = read.mapping_quality
		q_quality = read.query_alignment_qualities

		# Disregard any reads that don't have high mapping accuracy.
		if quality < 40:
			continue
		
		tmp = 0
		for i in range(tmp, len(positions)-1):

			# Check the probability that the base at this position is wrong.
			if q_quality[i] < 30 or sequence[i] == 'N':
				continue

			# Make sure that we compute just the specified region.
			if positions[i] >= end or positions[i] < start:
				break

			# Positions start at index 0 which is fine since reference also starts at 0.
			# Position number will be incremented for VCF file creation. 
			if positions[i] not in mutations:
				atcg = {'A': 0, 'T':0, 'C':0, 'G':0}
				atcg[sequence[i]] += 1
				mutations[positions[i]] = atcg
			else:				
				mutations[positions[i]][sequence[i]] += 1

	bamfile.close()	

	# No ref and separate mutation nucleotides fractions approach .
	#mutations = calculate_fractions(mutations, fastafile)
	
	# No ref and collected non-ref mutations fractions apprach.
	mutations, positions = calculate_fractions_overall(mutations, fastafile, chrom, re, rs)

	fastafile.close()

	return mutations, positions


def calculate_fractions(mutations, fastafile):
	"""
	Compares the mutations gathered from the BAM file to the reference FASTA file.
	Compiles the fractions of mutations (specific) vs overall number of reads.
	:param dict mutations: Dictionary of mutations with positions as key.
	:param file fastafile: Reference FASTA file object.
	:return: A newly compiled dictionary of mutations frequencies.
	"""

	# Calculate fractions and ignore reference nucleotide.
	for position, nucleotides in mutations.iteritems():

		# Open reference file for analysis.
		ref_n = fastafile.fetch("III", start=position, end=position+1)
		fractions = {'A': 0, 'T':0, 'C':0, 'G':0, 'N':0}

		# Get the total number of reads for specified position.
		all_reads = sum(nucleotides.values())

		for n in fractions.keys():
			if n is ref_n:
				continue
			else:
				fractions[n] = round(float(nucleotides[n])/all_reads, 4)

		mutations[position] = fractions

	return mutations


def calculate_fractions_overall(mutations, fastafile, chrom, re, rs):
	"""
	Compares the mutations gathered from the BAM file to the reference FASTA file.
	Compiles the fractions of mutations (not specific) vs overall number of reads.
	:param dict mutations: Dictionary of mutations with positions as key.
	:param file fastafile: Reference FASTA file object.
	:param float rs: Range start number.
	:param float re: Range end number.
	:return: A newly compiled dictionary of mutations frequencies.
	"""

	depth = 0
	positions = []

	# Create a VCF file of the positions with mutations so that it can be later used for VEP.
	# VCF will also be used for statistical analysis of the mutations present.
	f = open('viable_mutations_c_elegans.vcf', 'w')
	f.write('#CHROM 	POS 	ID 	REF 	ALT\n') 

	# Calculate fractions and ignore reference nucleotide.
	for position, nucleotides in mutations.iteritems():

		# Open reference file for analysis.
		ref_n = fastafile.fetch(chrom, start=position, end=position+1)
		muts = 0
		all_reads = sum(nucleotides.values())

		if all_reads < 100:
			mutations[position] = 0.0
			continue
		depth += all_reads

		for n in nucleotides.keys():
			if n is ref_n:
				continue
			else:
				muts += nucleotides[n]

		mutation_fraction = round(float(muts)/all_reads, 4)

		# Ignore any mutations that are over 90% because they might be just problems with reference. 
		if mutation_fraction > 0.9:
			mutations[position] = 0.0
			continue

		# Tell me the position of the mutation that had about 50% rate. 
		if mutation_fraction > rs/100 and mutation_fraction < re/100:
			positions.append(position)

			# Create a VCF file of the viable mutations.
			for k,v in mutations[position].items():
				if k is not ref_n:
					if float(v)/sum(mutations[position].values()) > 0.05:
						# VCF format for printing:
						# CHROM POS     ID        REF ALT 
						# read.get_reference_positions for some reason starts from 0 rather than 1. 
						f.write('%s %d C.E.0%d %s %s %f\n' %(chrom, position+1, position+1, ref_n, k, mutation_fraction*100))
		
		mutations[position] = mutation_fraction

	f.close()

	# Calculate the average depth of the read: How many data points did we get for each position?
	depth = float(depth) / len(mutations)
	print('The average depth of this segment reads was: %f' %depth)

	return mutations, positions


def vcf_analysis(vcf):
	"""
	Calculate the number of transversions and transitions in the vcf mutation file.
	Transitions (A<->G, C<->T)
	Transversions (A<->C, A<->T, C<->G, G<->T)
	:param string vcf: VCF file pathway.
	No return object.
	"""

	# Make a dictionary to accomodate all of the mutation types.
	mutations = {'A-G':0, 'C-T':0, 'A-C':0, 'A-T':0, 'C-G':0, 'G-T':0}

	with open(vcf, 'r') as f:
		# VCF format is [CHROM, POS, ID, REF, ALT] 
		for line in f:
			if not line.startswith('#'):
				ref = line.split()[3]
				alt = line.split()[4]
				key = ref + '-' + alt
				if key not in mutations:
					key = alt + '-' + ref
				mutations[key] += 1

	total = sum(mutations.values())
	if total:
		# Print out the analysis 
		print('Transitions    Occurence number    Percentage')
		print('A-G    %d     %f' %(mutations['A-G'], float(mutations['A-G'])/sum(mutations.values())))
		print('C-T    %d     %f' %(mutations['C-T'], float(mutations['C-T'])/sum(mutations.values())))
		print('\n\nTransversions  Occurence number    Percentage')
		print('A-C    %d     %f' %(mutations['A-C'], float(mutations['A-C'])/sum(mutations.values())))
		print('A-T    %d     %f' %(mutations['A-T'], float(mutations['A-T'])/sum(mutations.values())))
		print('C-G    %d     %f' %(mutations['C-G'], float(mutations['C-G'])/sum(mutations.values())))
		print('G-T    %d     %f' %(mutations['G-T'], float(mutations['G-T'])/sum(mutations.values())))
	else:
		print('No mutations detected.')

def test_files(args):
	"""
	Tests if the arguments specified are correct and can be used for the analysis.
	Exits the program if they are not.
	:param tuple args: Args object.
	No return values.
	"""

	# Check for bam file.
	try:
		bam = args.bamfile
	except:
		sys.exit('Error: Bam file not specified. Try: -bam bam_file_name_sorted.bam')

	# Check for reference file.
	try:
		ref = args.reference
	except:
		sys.exit('Error: Reference file not specified. Try: -ref reference.fasta')


	# Test if BAM file pathway is correct.
	try:
		pysam.AlignmentFile(args.bamfile, 'rb')
	except:
		sys.exit('Error: Unable to open BAM file.')

	# Test if reference file pathway is correct.
	try:
		pysam.FastaFile(args.reference)
	except:
		sys.exit('Error: Unable to open reference file.')


def test_numerical_args(args):
	"""
	Tests if the numerical arguments specified are correct and can be used for the analysis.
	Exits the program if they are not.
	:param tuple args: Args object.
	No return values.
	"""

	# Check for start position.
	try:
		start = int(args.start)
	except:
		sys.exit('Error: Start position not present or invalid.')

	# Check for end position.
	try:
		end = int(args.end)
	except:
		sys.exit('Error: End position not present or invalid.')

	# Make sure the start is smaller than the end position.
	if args.start >= args.end:
		sys.exit('Error: Invalid start and end positions specified.')

	# Test if the start and end positions make sense in context of the chromosome length.
	bam = pysam.AlignmentFile(args.bamfile, 'rb')
	header = bam.header
	for entry in (header['SQ']):
		if entry['SN'] == args.chromosome:
			length = entry['LN']

	if args.end > length:
		sys.exit('Error: End position exceeds the length of the chromosome.')

	# Test if ranges are in the appropriate range.
	try:
		rs = float(args.range_s)
	except:
		sys.exit('Error: Range start number invalid.')

	try:
		re = float(args.range_e)
	except:
		sys.exit('Error: Range end number invalid.')

	if args.range_s < 0 or args.range_e > 100:
		sys.exit('Error: Range numbers invalid. Provide range between 0 and 100 %.')


def main():

	# Parse the arguments. BAM, reference, start, and end positions need to be defined.
	parser = argparse.ArgumentParser()
	parser.add_argument('-bam', '--bamfile', help='BAM file path')
	parser.add_argument('-ref', '--reference', help='FASTA reference genome path')
	parser.add_argument('-s', '--start', help='Start position', type=int)
	parser.add_argument('-e', '--end', help='End position', type=int)
	parser.add_argument('-chr', '--chromosome', help='Chromosome number', type=str)
	parser.add_argument('-rs', '--range_s', help='Start of range of percentage of mutations.', 
						type=float, default=35.0)
	parser.add_argument('-re', '--range_e', help='End of range of percentage of mutations.', 
						type=float, default=65.0)

	args = parser.parse_args()

	# Make sure that user specified valid arguments. Close program if not.
	test_files(args)
	test_numerical_args(args)

	bam, ref, start, end, chrom = args.bamfile, args.reference, args.start, args.end, args.chromosome
	rs, re = args.range_s, args.range_e

	# Compile mutation frequencies. Return the frequencies at all positions and positions
	# 	of high frequency mutations.
	mutations, positions = find_saturation(bam, ref, start, end, chrom, rs, re)

	# Analyze the VCF file created from the saturation analysis.
	vcf_analysis('viable_mutations_c_elegans.vcf')


	overall_mutants = True

	# For the overall mutations approach.
	# Create a graph of reccurring mutations as each position.
	if overall_mutants:
		pd.DataFrame.from_dict(mutations,orient='index').plot()
		L=plt.legend()
		L.get_texts()[0].set_text('SNPs')
	else:
		# Show the plot for the separate nucleotides approach
		pd.DataFrame(mutations).T.plot()

	plt.ylim(top=1.0, bottom=0.0)
	step = (args.end-args.start)/10 # Show just 10 marks on the x axis.
	plt.xticks(np.arange(args.start, args.end, step=step))
	plt.show()
	# TODO: Save the map to a file
	# TODO: Go through the whole chromosome in chunks and save all of the plots in a 
	#		zipped folder. (Chunks of ~ 1,000,000 bp?)
	# TODO: Go through the whole chromosome in chunks and save the positions with mutations.
	#		Accumulate the dictionary of genes that have mutations and at the end rank them.


if __name__ == '__main__':
	main()