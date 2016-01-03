import argparse 
from contig_function import * 

parser = argparse.ArgumentParser(description="")

parser.add_argument('-i', '--input',
					dest = 'input',
					action = 'store',
					default = None, 
					help = 'Input linear file with filtered contigs information for each line')

parser.add_argument('-r', '--reference',
					dest = 'reference',
					action = 'store',
					default = None, 
					help = 'Input .fna file with the original contigs sequences.')

parser.add_argument('-o', '--output',
					dest = 'output',
					action = 'store',
					default = None,
					help = 'Output .fna file with the processed information')

options = parser.parse_args()

if __name__ == "__main__":

	if options.input is None:

		raise ValueError (">>> WARNING: Please, enter an input file to process. Write '-i (file)' <<<\n")


	else:
	
		infile = options.input

	if options.reference is None:

		raise ValueError (">>> WARNING: Please, enter a reference file to process. Write '-r (file)' <<<\n")

	else:

		reffile = options.reference

	if options.output is not None:

		fo = open(options.output, 'w')


	list_of_contigs = get_list_of_contigs(infile)

	for element in contig_fasta_sequence(reffile):

		if element[0] in list_of_contigs:

			if options.output is None:

				print ">"+element[0]+"\n"+element[1]+"\n"

			else:

				fo.write(">"+element[0]+"\n"+element[1]+"\n")


	if options.output is not None:

		fo.close()


