def contigs_info_generator(inputfile):

 	"""
 	A generator that computs a list of all the information of a file where the info is summaryzed
 	in a single line.

 	The list generated is the following:

 	Contig name - Num Genes - Total Genes Length (bp) - Num Exons - Total Exon Length (bp) - Mean Exon Length (bp) - Num Introns - Total Intron Length (bp) - Mean Intron Length (bp) - Intron Percent in Genes (w/ Introns) - Eukaryota Hits - Bacteria Hits - Archaea Hits - Virus Hits - Undefined Hits - Bad_hits )

 	"""


 	fh = open(inputfile, "r")

 	for line in fh:

 		if line.startswith("#"):

 			continue

 		else:

 			yield line.rstrip().lstrip().split("\t")

 	fh.close()


def contig_fasta_sequence(inputfile):

	'''
	A generator that yields the a list of two items, the first one is the name of the contig and the second one is the nucleotid sequence. It needs as an input a FASTA file of the masked contigs. 

	'''

	fh = open(inputfile, 'r')

	number_of_lines = 0 

	for line in fh:

		number_of_lines += 1

		if line.startswith(">"):

			if number_of_lines > 1:

				yield [contig_name, sequence]

			contig_name = line.lstrip(">").rstrip()

		else:

			sequence = line.rstrip().lstrip()

	yield [contig_name, sequence]

	fh.close()


def get_contig_name(input_list):

	'''
	Return the name of the contig as a string.

	'''

	return input_list[0]

def get_contig_len(input_list):

	'''
	Return the length of the contig as an integer.

	'''
	
	return int(get_contig_name(input_list).split("_")[3])

def get_number_genes(input_list):

	'''
	Return the number of genes that a contig has as an integer.

	'''

	return int(input_list[1])

def get_total_genes_len(input_list):

	'''
	Return the total length of genes within the contig as an integer in pb.

	'''

	return int(input_list[2])

def get_number_exons(input_list):

	'''
	Return the number of exons that a contig has as an integer.

	'''

	return int(input_list[3])

def get_total_exon_len(input_list):

	'''
	Return the total length of exons within the contig as an integer in pb.

	'''

	return int(input_list[4])

def get_mean_exon_len(input_list):

	'''
	Return the mean length of exons within the contig as an integer in pb.

	'''

	return float(input_list[5])

def get_number_intron(input_list):

	'''
	Return the number of introns that a contig has as an integer.

	'''

	return int(input_list[6])

def get_total_intron_len(input_list):

	'''
	Return the total length of introns within the contig as an integer in pb.

	'''

	return int(input_list[7])

def get_mean_intron_len(input_list):

	'''
	Return the mean length of introns within the contig as an integer in pb.

	'''

	return float(input_list[8])

def get_intron_percent_of_genes(input_list):

	'''
	Return the introns' one percent respect the genes within the contig that have introns.


	'''

	return float(input_list[9].rstrip("%"))/100

def get_euk_hits(input_list):

	'''
	Return an integer of the number of eukaryotic hits. 

	'''

	return int(input_list[10])

def get_bact_hits(input_list):

	'''
	Return an integer of the number of bacterial hits. 

	'''

	return int(input_list[11])

def get_virus_hits(input_list):

	'''
	Return an integer of the number of viral hits. 

	'''

	return int(input_list[12])

def get_arch_hits(input_list):

	'''
	Return an integer of the number of archaea hits. 

	'''

	return int(input_list[13])

def get_undef_hits(input_list):

	'''
	Return an integer of the number of undefined hits. 

	'''

	return int(input_list[14])

def get_bad_hits(input_list):

	'''
	Return an integer of the number of bad hits. 

	'''

	return int(input_list[15])


def intron_exon_coverage(input_list):

	''' 
	Returns a two items' tuple where the first item is the intron coverage and the second one
	is the exon coverage.

	'''

	intron_cov = float(get_total_intron_len(input_list))/get_contig_len(input_list)

	exon_len = float(get_total_exon_len(input_list))/get_contig_len(input_list)

	return (intron_cov, exon_len)

def get_list_of_contigs(inputfile):

	'''
	From an inputfile with one line per contig and it's information, this function returns a list of all the contigs present in the file.

	'''

	contigs_list = []

	for contig in contigs_info_generator(inputfile):

		contigs_list.append(get_contig_name(contig))

	return contigs_list
	


