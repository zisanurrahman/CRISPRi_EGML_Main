import argparse
import re 
import numpy as np
import csv

# Definitions and constants --------------------------------------------------------------------
min_uniqueness = 12
start_offset = -300			# Upstream of ATG, change this number as you like
end_offset = 120			# Downstream of ATG, change this number as you like
path = "k52.gb"	# Name of the input file
max_repetition = 4

parser = argparse.ArgumentParser(description='Find interesting sequences in a .gb file')
parser.add_argument('path', metavar='genebank_file', type=str, help='Input genebank file.')
parser.add_argument('start_offset', metavar='so', type=int, help='Upstream offset.')
parser.add_argument('end_offset', metavar='eo', type=int, help='Downstream offset.')
parser.add_argument('max_repetition', metavar='mr', type=int, help='Max number of consecutive repeated base pairs.')

try:
	args = parser.parse_args()
	path = args.path
	start_offset = args.start_offset * -1
	end_offset = args.end_offset
	max_repetition = args.max_repetition
except:
	print("Some arguments not provided. Using defaults.")

print ("MAX REPETITION: ", max_repetition)

end_sequence = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC'.lower()
end_83 = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT'.lower()

nucleotides = ['a','c','t','g']

complement = {
	'a': 't',
	'c': 'g',
	't': 'a',
	'g': 'c'
}

# Helper methods --------------------------------------------------------------------------------

'''
	reverse_complement

	Input: A sequence of nucleotides
	Output: The reverse complement of the sequence

	Description: Uses a dictionary lookup to replace every letter (nucleotide) with its complement,
	then returns the reverse of that string.
'''
def reverse_complement(s):
	return ''.join([complement[c] for c in s][::-1])


'''
	parse_annotation

	Input: A single line of text from the .gb file which contains a gene annotation.
	Output: A dictionary containing:
		-The start position of the gene
		-The end position of the gene
		-A boolean, "complement", which is true if the gene is on the opposite strand

	Description: Checks if the line contains the word 'complement', which indicates that the gene
	is found in the opposite strand. Then, parses the line to find the start and end of the gene.
	If the line given was not actually an annotation, then split() will not return the correct
	number of "parts", and we return None to indicate that there was no gene.
'''
def parse_annotation(line):
	if 'complement' in line:
		complement = True
		bp_range = line[11:-1]	# This range turns complement(start..end) into start..end
	else:
		complement = False
		bp_range = line

	bp_range = bp_range.split('..')
	if len(bp_range) == 2:
		start = int(re.sub(r'[<>]', '', bp_range[0]))	# We remove the '>' or '<' which indicate partial features
		end = int(re.sub(r'[<>]', '', bp_range[1]))
		return {
			'start': start,
			'end': end,
			'complement': complement
		}
	else:
		return None


'''
	parse_file

	Input: A string, the path to the file we wish to parse.
	Output: The entire genome as a string, with one character representing each nucleotide,
		The annotations as a list of dictionaries, with one dictionary per gene.

	Description: Opens the file to read. Once we see "ORIGIN" we know we are in the genome and
	can directly append each character read to our genome string. Otherwise, we check if the
	line we are reading contains the word 'gene', which generally indicates an annotation in the
	.gb file format. If so, we pass the line to parse_annotation to parse the information about
	the gene annotation.
'''
def parse_file(filename):
	genome = ''
	annotations = []
	with open(filename) as f:
		in_genome = False
		for line in f:
			line = line.strip()
			if line == "ORIGIN":
				in_genome = True
			elif in_genome:
				genome = genome + ''.join([c for c in line if c in nucleotides])
			else:
				line = line.split()
				if line[0] == 'gene':
					annotation = parse_annotation(line[1])
					if annotation is not None:
						annotations.append(annotation)
								
	return genome, annotations


'''
	gen_position_vector

	Input: The genome string, and a list of all the gene annotations.
	Output: A vector for each strand, each vector containing a value for every
	nucleotide in the genome, indicating whether or not that position is part of a gene.

	Description: First, creates an array of 0s with the same length as the genome, with one
	array per strand. Then, goes through every annotation and flips positions "on" if they are
	within the specified range of the gene's start. Each annotation only affects one vector,
	determined by whether or not the annotations has the "complement" flag.

	Then, for every gene we go back and mark positions "off" again if they are between end_offset
	nucleotides after the gene's start, and the end of the gene. Again, this works only on one
	strand's vector.
'''
def gen_position_vector(genome, annotations):
	_min = 0
	_max = len(genome)

	print("Creating empty position vectors...")
	vec53 = np.zeros(_max).astype(int)
	vec35 = np.zeros(_max).astype(int)
	vec53comp = np.ones(_max).astype(int)
	vec35comp = np.ones(_max).astype(int)

	print("Parsing gene range annotations...")
	for a in annotations:
		if a['complement']:
			start = min(_max, a['end'] - start_offset)
			end = max(_min, a['end'] - end_offset)
			vec35[end:start] = 1
			
			# We don't want any positions too deep 'inside' the gene to be considered
			vec35comp[a['start']:end] = 0

			if end > start:
				print("oh no!")
		else:
			start = max(_min, a['start'] + start_offset)
			end = min(_max, a['start'] + end_offset)
			vec53[start:end] = 1
			
			# We don't want any positions too deep 'inside' the gene to be considered
			vec53comp[end:a['end']] = 0

			if end < start:
				print("oh no!")

	print("Finalizing gene range vectors...")
	vec53 = np.logical_and(vec53, vec53comp, dtype=int)
	vec35 = np.logical_and(vec35, vec35comp, dtype=int)

	print(sum(vec53) + sum(vec35), " active positions.")

	return vec53, vec35


'''
	get_bps_in_range

	Input: The genome, and a vector indicating which positions in the genome we want to
	search for matches in.

	Output: The slices of the genome in which we want to look for matches, as separate
	lists of nucleotides and their positions. 

	Description:
	This results in a new smaller "genome" with no visible gaps. Thus, when using this smaller 
	genome, you must check if you crossed a "gap" by checking the positions of the nucleotides 
	you're using.

	For example, imagine the input: (CCGCGG, 110011). This would output (CCGG, [1,2,5,6]).
	If you were looking for the pattern "CG", you would find it, but would notice that the
	C was position 2, and G was position 5, and thus the match is "invalid" since there would
	be other nucleotides between those nucleotides in the real genome.

	The reason for this method is that the genome itself is huge and expensive to look through;
	it is much faster to find all matches in this smaller genome, and then filter out the false
	positives that result due to jumping across different parts of the genome.
'''
def get_bps_in_range(genome, vec):
	print("Creating sub-genome based on ranges...")
	nucleotides = ''
	positions = []
	for i in range(len(genome)):
		if vec[i]:
			nucleotides += genome[i]
			positions.append(i)

	return nucleotides, positions
			

'''
	get_candidates

	Input: A list of nucleotides, and which position in the genome these nucleotides
	are.

	Output: A list of dictionaries for every match. Each dictionary contains the 20bp match,
	as well as the start position of this nucleotide string.

	Description: A match is found by searching for the regular expression of: 
		cc(any nucleotide)(20 more nucleotides)
	
	Before accepting a match we must make sure that there are no "gaps" in the positions, since
	we are working with a reduced subset of the original genome. This is done by making sure that
	the start position and end position have a difference of 23 nucleotides (CCX + 20).
	
'''

def get_candidates(nucleotides, positions):
	print("Searching for 20bp candidates...")
	pattern = re.compile("cc\w\w{20}")

	results = []

	for m in pattern.finditer(nucleotides):
		start = m.start()
		# Make sure we don't keep false positives which result from "jumping" across active range borders
		# For example, if [1..15] are marked active and [17..30] are, we don't want to keep matchines
		# which rely on having a nucleotide from both positions 15 and 17.
		if (positions[start + 23] - positions[start]) == 23:
			nucleotides = m.group()[3:]

			match = {
				'start': positions[start] + 3,
				'nucleotides': nucleotides,
				'strand': ''
			}

			results.append(match)

	print(len(results), " candidates.")
	return results


'''
	filter_duplicates

	Input: The whole genome sequence (seq), and a list of candidate sequences.

	Output: A list of unique candidate sequences.

	Description: A candidate sequence is considered unique if:
		-The entire sequence appears exactly once in the genome
		-The first [19, 18 ... min_uniqueness] nucleotides appear as a sequence 
		 exactly once in the genome.

'''

def is_repetitive(seq):
	# 1 is subtracted from max_repetition to make the regex capture the correct number
	regex = r'(\w)\1{' + str(max_repetition - 1) + r'}'
	return(len(re.findall(regex, seq)) > 0)

def filter_duplicates(seq, candidates, reverse=False):
	uniques = []
	index = 1
	total = len(candidates)

	print("Checking uniqueness of candidate nucleotide sequences.")
	for candidate in candidates:
		nucleotides = candidate['nucleotides']

		if not is_repetitive(nucleotides):
			if reverse:
				nucleotides = reverse_complement(nucleotides)

			if index % 1000 == 0:
		   	 print('Checking candidate sequence {} of {}.'.format(index, total))

			index += 1

			# For a given string, this will check if all chars are unique,
			# then the first 19, then the first 18... down to min_uniqueness
			unique = True
			for i in range(20 - min_uniqueness):
				if seq.count(nucleotides[:min_uniqueness+i]) != 1:
					unique = False
					break
			if unique:
				uniques.append(candidate)
		else:
			print("skipping repetitive sequence: ", nucleotides)


	for unique in uniques:
		if reverse:
			unique['strand'] = '-' 
		else:
			unique['strand'] = '+'

	print("After filtering by uniqueness: ", len(uniques))	
	return uniques


def reverse_comp_list(seqs):
    return[{'start': seq['start'], 'nucleotides': reverse_complement(seq['nucleotides']), 'strand': seq['strand']} for seq in seqs]


def add_35bp(seqs):
    return[{'start': seq['start'], 'nucleotides': seq['nucleotides'] + end_sequence, 'strand': seq['strand']} for seq in seqs]
	
def add_83bp(seqs):
    return[{'start': seq['start'], 'nucleotides': seq['nucleotides'] + end_83, 'strand': seq['strand']} for seq in seqs]

# Main program ----------------------------------------------------------------------------------


print("Loading genome and annotations...")
genome, annotations = parse_file(path)

vec53, vec35 = gen_position_vector(genome, annotations)

nucleotides, positions = get_bps_in_range(genome, vec53)
nucleotides35, positions35 = get_bps_in_range(genome, vec35)

candidates = get_candidates(nucleotides, positions)
candidates35 = get_candidates(nucleotides35, positions35)

uniques = filter_duplicates(genome, candidates)
uniques35 = filter_duplicates(genome, candidates35, reverse=True)

uniques = reverse_comp_list(uniques)
uniques = uniques + uniques35

# save as csv
print("Saving unique results to csv")
with open("results.csv", 'w') as f:
	writer = csv.DictWriter(f, ['start', 'nucleotides', 'strand'])
	writer.writeheader()
	for data in add_35bp(uniques):
		writer.writerow(data)

# save fasta
with open("results.fasta", 'w') as f:
	for data in add_83bp(uniques):
		line1 = '>' + str(data['start']) +  str(data['strand'])
		line2 = data['nucleotides']
		f.write(line1 + '\n')
		f.write(line2 + '\n')

print("Saved results.csv and results.fasta")

