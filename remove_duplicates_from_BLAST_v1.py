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

options = parser.parse_args()

###### Opening the file ######

if options.input is None:

	raise ValueError (">>> WARNING: Please, enter an input file to process. Write '-i (file)' <<<\n")

else:
	
	infile = options.input

fh = open(infile, "r")

blast_entries_dict = {}

list_of_unique_entries = []

h = 0

for line in fh:

	tmp_array = line.rstrip().split("\t")

	if tmp_array[0] not in list_of_unique_entries: 

		if h > 0:

			blast_entries_dict[entry_name] = info_list

		list_of_unique_entries.append(tmp_array[0])

		info_list = []

	entry_name = tmp_array[0]

	info_list.append(tmp_array[1:])

	h += 1

blast_entries_dict[entry_name] = info_list

fh.close()

if options.output is None:

	fh2 = open(infile.rstrip('.txt')+"_non_duplicates.txt", 'w')

else:

	fh2 = open(options.output, "w")

for entry in list_of_unique_entries:

	if entry in blast_entries_dict.keys() and len(blast_entries_dict[entry]) > 1:

		blast_entries_dict[entry].sort(key=lambda x: float(x[2])/float(x[11]), reverse=True)

	fh2.write(entry+"\t"+"\t".join(blast_entries_dict[entry][0])+"\n")

fh2.close()








