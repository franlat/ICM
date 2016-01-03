import argparse 
from numpy import array
from numpy import *

parser = argparse.ArgumentParser(description="")

parser.add_argument('-i', '--input',
					dest = 'input',
					action = 'store',
					default = None, 
					help = 'Input linear file with an exon per line from the AUGUSTUS gene prediction')

parser.add_argument('-o', '--output',
					dest = 'output',
					action = 'store',
					default = None,
					help = 'Output linear file with the processed information')

parser.add_argument('-b', '--blast',
					dest = 'blast', 
					action = 'store', 
					default = None,
					help = 'File with the taxonomy of each gene')

parser.add_argument('-id', '--identity',
					dest = 'identity',
					action = 'store',
					type = float,
					default = 70,
					help = 'Minimum percent of identity required in the blast file to be considered correct')

parser.add_argument('-co', '--coverage',
					dest = 'coverage',
					action = 'store',
					type = float,
					default = 60,
					help = 'Minimum percent of coverage required in the blast file to be considered correct')

options = parser.parse_args()

###### Opening the file ######

if options.input is None:

	raise ValueError (">>> WARNING: Please, enter an input file to process. Write '-i (file)' <<<\n")

elif options.input[-3:] != "gff":
	
	raise ValueError (">>> WARNING: Please, enter a .gff file to process <<<")

else:
	
	infile = options.input

fh = open(infile, "r")

###### Defining the counters ######

contigs = 0 
genes_count = 0 
introns = 0
introns_len = 0 
exons_len = 0 
exons = 0 

###### Defining the strings ######

protein_seq = ''

###### Defining the lists ######

genes = []
final_genes = []

###### Definin the final dictionary ######

genes_info = {}
final_dict = {}

###### Processing ######

for line in fh:

	if line.startswith("# -"):

		if contigs > 0:

			final_genes.append(contig_name + genes)

		start = line.index("NODE")
		end = line.index(")")
		contig_name = [line[start:end]]
		contigs += 1
		genes = []

	elif line.startswith("# start"):

		if genes_count > 0:

			if introns == 0 :

				introns_mean = float(0)
				introns_percent = float(0)

			else:

				introns_mean = float(introns_len)/introns
				introns_percent = float(introns_len)/gene_len

			genes_info[gene_name] = [gene_len, exons, exons_len, float(exons_len)/exons, introns, introns_len, introns_mean, introns_percent, protein_seq]

		gene_name = line.split()[3]	
		genes.append(gene_name)

		# Reset the counters for the next gene

		introns = 0
		introns_len = 0 
		exons = 0 
		exons_len = 0  
		protein_seq = ''
		genes_count += 1


	elif line.startswith("# namgene"):

		genes.append("No predicted genes")

	elif line.startswith("NODE"):

		if line.split()[2] == "gene":

			gene_len = int(line.split()[4]) - int(line.split()[3])

		elif line.find("intron") != -1:

			introns += 1
			introns_len += int(line.split()[4]) - int(line.split()[3])

		elif line.find("CDS") != -1:

			exons += 1
			exons_len += int(line.split()[4]) - int(line.split()[3])

	elif line.startswith("# protein"):

		start_aa = line.index("[") + 1
		protein_seq += line.rstrip().rstrip("]")[start_aa:]

	elif line.startswith("# ") and line.isupper():

		if line.rstrip().endswith("]") != True:

			protein_seq += line.rstrip()[2:]

		else:
			
			protein_seq += line.rstrip()[2:-1]

	else:

		continue


final_genes.append(contig_name + genes)

if introns == 0:

		introns_mean = 0
		introns_percent = 0 

else:

		introns_mean = float(introns_len)/introns
		introns_percent = float(introns_len)/gene_len

genes_info[gene_name] = [gene_len, exons, exons_len, float(exons_len)/exons, introns, introns_len, introns_mean, introns_percent, protein_seq]

fh.close()




###### Opening the BLAST file ######


if options.blast is None:

	raise ValueError (">>> WARNING: Please, enter an input file to process. Write '-b (file)' <<<\n")

else:

	blast_file = options.blast

fh3 = open(blast_file, "r")


###### Defining the counters ######

Eukaryota = 0 
Bacteria = 0
Archaea = 0
Virus = 0 
Undefined = 0
Bad_hits = 0 

h = 0



###### Defining the list and dict of genes ######


genes_bl = []  # A list to take into acount the genes with more than one exon

tmp_bl = []

genes_class = {}

final_genes_bl = {} # The final list to display the results

###### How many exons a gene has and where do they belong ######

for line in fh3:

	tmp_list = line.split("\t") # Split the line by columns
	index = tmp_list[0].index(".")  # We llook for the index in where the gene name ends

	if tmp_list[0][:index] not in genes_bl:  #first time we see the gene in the list

		if h > 0:   # to avoid adding the first line when there is nothing computed yet

			final_genes_bl[gene_name] = [Good_hits, Bad_hits]		

		genes_bl.append(tmp_list[0][:index])   # We append the gene name into a the list of 'first time'

		Good_hits = 0   # We reset all the variables to 0
		Bad_hits = 0 
	
	gene_name = tmp_list[0][:index]  # We add the name to a new variable to be able to identify it later before it is rewritten.

	if tmp_list[2] >= options.identity and float(tmp_list[3])*100/int(tmp_list[12]) >= options.coverage:

		Good_hits += 1

	else:

		Bad_hits += 1

	h += 1

final_genes_bl[gene_name] = [Good_hits, Bad_hits] #since no new gene will be found at the end of the file, the last gene wouldn't appear in the results

fh3.close()


###### Computing the contigs' data ######

contigs_info = {}

for element in final_genes:

	if element[1] != 'No predicted genes':

		tmp_hits_list = []
		tmp_gene_info = []

		for i in range(1, len(element)):   # to sum up the hits of a contig and the info

			if element[i] in genes_info.keys():

				tmp_info = genes_info[element[i]]

				# genes_info[gene_name] = [gene_len, exons, exons_len, float(exons_len)/exons, introns, introns_len, introns_mean, introns_percent, protein_seq]
				
				if tmp_info[4] == 0:

					tmp_gene_info.append([tmp_info[0], tmp_info[1], tmp_info[2], tmp_info[4], tmp_info[5], 0])

				else:

					tmp_gene_info.append([tmp_info[0], tmp_info[1], tmp_info[2], tmp_info[4], tmp_info[5], tmp_info[0]])

			else:

				tmp_gene_info.append([0]*6)

			if element[i] in final_genes_bl.keys():	

				tmp_hits_list.append(final_genes_bl[element[i]])

			else:

				tmp_hits_list.append([0]*2)


		first_list = sum(array(tmp_gene_info), 0)

		second_list = sum(array(tmp_hits_list), 0)

		no_genes = len(element[1:])

		if first_list[3] == 0:

			intr_mean = float(0)
			intr_percent = float(0)

		else:

			intr_mean = float(first_list[4])/first_list[3]
			intr_percent = float(first_list[4])/first_list[5]

		# Number of genes - length of genes - number of exons - total exon length - mean exon length - number of introns - total intron length - mean introng length - % of introns in genes (with introns) - HITS

		contigs_info[element[0]] = append([no_genes, first_list[0], first_list[1], first_list[2], float(first_list[2])/no_genes, first_list[3], first_list[4], intr_mean, intr_percent], second_list).tolist()


###### Pritning the results ######

if options.output is None:

	fh2 = open("contigs_parser_UNIGENE_results_v1.txt", 'w')

else:

	fh2 = open(options.output, "w")


fh2.write("This file is generated from a .gff file of the AUGUSTUS prediction and the results of a BLAST process.\n")
fh2.write("Lines starting with # are the name of the contig followed by a serie of data:\n")
fh2.write("\n# Contig name - Num Genes - Total Genes Length (bp) - Num Exons - Total Exon Length (bp) - Mean Exon Length (bp) - Num Introns - Total Intron Length (bp) - Mean Intron Length (bp) - Intron Percent in Genes (w/ Introns) - Good Hits - Bad_hits )\n")
fh2.write("Lines starting with > are the genes followed by a serie of data:\n")
fh2.write("\n>GeneName - Gene Length (bp) - Num Exons - Total Exon Length (bp) - Mean Exon Length (bp) - Num Introns - Total Intron Length (bp) - Mean Intron Length (bp) - Intron Percent in Gene - Good Hits - Bad_hits\n")
fh2.write(">>> (Bad Hits the ones that do not match the required Identity and Coverage) <<<\n")
fh2.write("\nThe following line is the Protein Sequence.\n\n")


for element in final_genes:

	if element[0] in contigs_info.keys():

		contig_header = "\t{0[0]:.0f}\t{0[1]:.0f}\t{0[2]:.0f}\t{0[3]:.0f}\t{0[4]:.2f}\t{0[5]:.0f}\t{0[6]:.0f}\t{0[7]:.2f}\t{0[8]:.2%}\t{0[9]:.0f}\t{0[10]:.0f}"

		fh2.write("\n# "+element[0]+contig_header.format(contigs_info[element[0]])+"\n\n")

	else:

		contig_header = "\t{0[0]:.0f}\t{0[1]:.0f}\t{0[2]:.0f}\t{0[3]:.0f}\t{0[4]:.2f}\t{0[5]:.0f}\t{0[6]:.0f}\t{0[7]:.2f}\t{0[8]:.2%}\t{0[9]:.0f}\t{0[10]:.0f}"

		fh2.write("\n# "+element[0]+contig_header.format([0]*11)+"\n\n")
		
	if element[1] == "No predicted genes":

		fh2.write(element[1]+"\n")

	else:

		for i in range(1, len(element)):

			if element[i] in genes_info.keys():
				
				if element[i] in genes_bl:

					fh2.write(">"+element[i]+"\t{0[0]}\t{0[1]}\t{0[2]}\t{0[3]:.2f}\t{0[4]}\t{0[5]}\t{0[6]:.2f}\t{0[7]:.2%}\t{1[0]}\t{1[1]}\n".format(genes_info[element[i]], final_genes_bl[element[i]]))

				else:

					fh2.write(">"+element[i]+"\t{0[0]}\t{0[1]}\t{0[2]}\t{0[3]:.2f}\t{0[4]}\t{0[5]}\t{0[6]:.2f}\t{0[7]:.2%}\t{1[0]}\t{1[1]}\n".format(genes_info[element[i]], [0]*2))

				fh2.write(genes_info[element[i]][8]+"\n\n")

			else:

				fh2.write(">"+element[i]+"\n")

			# if element[i] in final_genes_ex.keys():

				# fh2.write("- ".rjust(5)+element[i].ljust(6)+blast_info.format(final_genes_ex[element[i]]))

			# else:

				# fh2.write("- ".rjust(5)+element[i].ljust(6)+blast_info.format([0]*6))

fh2.close()















