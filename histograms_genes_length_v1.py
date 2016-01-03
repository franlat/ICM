from contig_function import * 
import argparse
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math

parser = argparse.ArgumentParser(description="")

parser.add_argument('-i', '--input',
					dest = 'input',
					action = 'store',
					default = None, 
					help = 'Input linear file with only the contigs data from the ContigsParser')

parser.add_argument('-o', '--output',
					dest = 'output',
					action = 'store',
					default = None,
					help = 'Output linear file with the processed information')

options = parser.parse_args()

if __name__ == "__main__":

	if options.input is None:

		raise ValueError (">>> WARNING: Please, enter an input file to process. Write '-i (file)' <<<\n")

	else:
	
		infile = options.input


	total_num_genes_list = []
	total_len_contigs_list = []

	for element in contigs_info_generator(infile):

		num_genes = get_number_genes(element)
		contig_len = get_contig_len(element)

		total_num_genes_list.append(num_genes)
		total_len_contigs_list.append(contig_len)

	if options.output is None:



		plt.hist(total_num_genes_list, bins=max(total_num_genes_list), histtype='step')
		plt.xlabel("Number of genes")
		plt.ylabel("Number of Contigs")
		plt.title('Histogram')
		plt.grid(True)
		plt.show()

		plt.hist(total_len_contigs_list, bins=int(math.floor(max(total_len_contigs_list)/1000))+1, histtype='step')
		plt.xlabel("Contig Length (bp)")
		plt.ylabel("Number of Contigs")
		plt.title('Histogram')
		plt.grid(True)
		plt.show()

	else:

		fh2 = open(options.output, 'w')

		fh2.write("Number Genes\tContig Len\n")

		for i in range(0, len(total_num_genes_list)-1):

			fh2.write(str(total_num_genes_list[i])+"\t"+str(total_len_contigs_list[i])+"\n")

		fh2.close()
