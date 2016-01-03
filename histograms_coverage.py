from contig_function import * 
import argparse
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

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


	total_exon_coverage_list = []
	total_intron_coverage_list = []

	for element in contigs_info_generator(infile):

		in_ex_list = intron_exon_coverage(element)

		total_exon_coverage_list.append(in_ex_list[1]*100)
		total_intron_coverage_list.append(in_ex_list[0]*100)

	if options.output is None:

		plt.hist(total_intron_coverage_list, bins=500, histtype='step')
		plt.xlabel("Introns' Coverage %")
		plt.ylabel("Number of Contigs")
		plt.title('Histogram')
		plt.grid(True)
		plt.show()

		plt.hist(total_exon_coverage_list, bins=500, histtype='step')
		plt.xlabel("Exons' Coverage %")
		plt.ylabel("Number of Contigs")
		plt.title('Histogram')
		plt.grid(True)
		plt.show()

	else:

		fh2 = open(options.output, 'w')

		fh2.write("Intron Coverage\tExon Coverage\n")

		for i in range(0, len(total_intron_coverage_list)-1):

			fh2.write(str(total_intron_coverage_list[i])+"\t"+str(total_exon_coverage_list[i])+"\n")

		fh2.close()









