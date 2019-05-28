# UCSB
README for mutation_saturation.py 
Owner: Ada Madejska, UCSB, J.Rothman's lab 2019

Description:
	This program can be used for looking for SNPs frequency across multiple organisms
	of the same species. It creates a graph of specified region from the genome and 
	looks at the frequency at which the nucleotides differ from the reference at each
	position. Lastly, it creates a VCF file with all viable mutations. 

Input:
	Sorted and indexed BAM file (1) with all sequenced and aligned genomes for study
	Sorted and indexed FASTA reference file with which the BAM files were aligned.
	Start and end position (soon to be discontinued and replaced by the chromosome number).

Output:
	A graph showing the fractions of mutations happening at each position in the 
	specified regions.
	List of genes affected by the most frequent SNPs with number of mutations.
	A file with viable mutations in a VCF (Variant Calling Format) format. 

How to use it:
	python mutation_saturations.py -bam sorted_bam.bam -ref sorted_ref.fa -s start_pos
	-e end_pos -chr chromosome_number

How it works:
	mutation_saturation.py takes the sorted and indexed BAM file and goes through 
	the specified region of the chromosome (or the whole chromosome later in the 
	development). It counts all of the datapoints it gets from all of the reads 
	for each chromosome position (creates a dictionary where keys are chrom positions
	and values are the bases within another dictionary).
	After it is done with the the creation of the data dictionary, it reads the 
	reference genome to see which nucleotide at the specified position is 
	the reference nucleotide. It sums the other nucleotides together and 
	calculates a fraction (mutations/total reads). It the stores that information
	in the dictionary. 
	In addition to calculating fraction, the program also keeps track if the 
	fraction is between 35% and 60%. If the boundaries are hit,
	we say that the mutation is important and we store the position of it in a 
	list.
	Next, it creates a VCF (variant calling format) file with the list of all 
	viable mutations. It provides the chromosome number, ID, position of the mutation,
	reference nucleotide, mutated nucleotide, percentage of mutations vs all reads
	in this particular position. 
	The program also outputs the graph showing the specified chromosome region
	and the mutation frequency at each nucleotide position. It saves the VCF under 
	user-specified name. It also outputs statistical analysis of transitions vs 
	transversions within the specified region (See notes). 
	Later, a Variant Effect Predictor (VEP) functionality will be added. 
	As for May 2019, to collect data on mutations that result in important protein 
	coding changes in the genome, the user needs to manually run VEP on the 
	program-generated VCF file. 

Notes:
	Multiple concerns about the logic and execution of the program has been discussed
	and resolved. 
	1. Bell curved shaped mutation frequencies.
		Without accounting for the mapping quality of the reads and overall quality of the alignment,
		the graphs of mutations frequencies has shown a very distinct bell-shaped curves.
		This problem has been resolved by stringent quality control of the data.
		Many reads are being discarded if their mapping quality is low (less than 40).
		Many basepair alignment are being discarded if their individual quality score is low. 
		Additionally, we are not counting any 'N' nucleotides in the mutation frequency.
	2. The amount of reads per position in the chromosome after the quality control.
		After we discard a lot of low quality nucleotides and reads, how much is left? 
		Is there enough to make the mutation data believable? If the total number of data points
		for the particular position is only 1 or 2, a simple sequencing error might show 100% 
		mutation frequency at this position.
		To address this problem, we are not concidering any mutation frequencies that 
		have less than 100 datapoints. If there's less than 100 points at a particular position,
		this position defaults to 0% mutation frequency. 
		Additionally, the program outputs into the console the average amount of data points
		for each position so the user can decide if they believe the analysis. 
	3. Transitions vs transversions. Why is it importnat?
		You can divide nucleotides ATCG into two big groups : purines (structure has two rings) 
		(A and G) and pyrimidines (structure has only one ring) (T and C). Transitions are interchanges
		of the two ring structure into a one ring structure and vice versa. Transversions are 
		changes of nucleotides into the same group (two ring into two ring structure or one ring
		into one ring structure). Because of their structure, transitions are less likely to result 
		in amino acid substitutions (due to "wobble" - permition of  several types of non-Watson–Crick
		 base pairing to occur at the third codon position), and are therefore more likely to persist as 
		"silent substitutions" in populations as single nucleotide polymorphisms (SNPs). They are usually
		generated more than the structural change transversions. 
		The generation of transversions during replication requires much greater distortion of
		the double helix than does the production of transition mutations.
