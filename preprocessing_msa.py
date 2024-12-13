import os
import sys
import itertools
import collections
#from Bio import SeqIO
import pdb

# make amino acids set 
aa = "ARNDBCEQZGHILKMFPSTWYV"
amino_acids = set()
for c in aa:
	amino_acids.add(c)

# I need the named tuple to count amino acids at each position 
# and order them by the count easily
AA = collections.namedtuple('AA', 'name count')


def read_fasta(fasta_file_path):
    """
    read fasta file and return a dict of all sequences

    :param fasta_file_path: path to fasta file
    :return: a dictionary of seq_name:seq
    """

    sequences = dict()

    if not os.path.exists(fasta_file_path):
        logging.error("file {} does not exist".format(fasta_file_path))
        sys.exit()

    if fasta_file_path.endswith("gz"):
        fasta_file = gzip.open(fasta_file_path, "rt")
    else:
        fasta_file = open(fasta_file_path, "r")

    seqs = []
    seq_name = ""
    for line in fasta_file:
        line = line.strip()
        if not line:  # empty line
            continue

        if line.startswith(">"):
            if len(seqs) != 0:  # there was a sequence before
                sequences[seq_name] = "".join(seqs)
                seq_name = line[1:]
                seqs = []
            else:
                seq_name = line[1:]
        else:
            seqs.append(line)

    if seqs:
        sequences[seq_name] = "".join(seqs)

    return sequences


def return_pos_counts(reads, percentage, dna=False):
	'''
	This function takes the aligned reads dictionary and returns a list of lists, each list for a position in the aligned proteins
	and the distribution of the amino acids at that position in the form of AA namedtuples

	:params reads: dictionary of sample IDs as keys and protein sequence as values
	'''
	position_counts_dict = []
	position_counts_dict.append(dict())
	len_reads = len(reads)

	for i in range(len(reads[list(reads.keys())[1]])):
		position_counts_dict.append(dict())
		for read in reads.values():
			if read[i] not in position_counts_dict[i]:
				position_counts_dict[i][read[i]] = 1
			else:
				position_counts_dict[i][read[i]] += 1

	position_counts = []
	for pos, count_list in enumerate(position_counts_dict):
		position_counts.append([])  # list of AA objects
		for aa in count_list.keys():
			if dna:
				# skipping samples that have an ambiguos letter at that position
				if aa in {"-", "*", 'B', 'D','H', 'K', 'M', 'N', 'R', 'S','W', 'Y'}:
					continue
			else:
				# skipping samples that have a - or * at that position
				if aa in {"-", "*"}:
					continue

			# discard letters with counts less than percentage provided
			if count_list[aa] < percentage*len_reads:
				continue
			else:
				position_counts[pos].append(AA(name=aa, count=count_list[aa]))

	for count_list in position_counts:
		count_list.sort(key = lambda x: x.count, reverse = True)

	return position_counts


def generate_vcf_files(possible_pairs, position_counts, reads):
	'''
	Generates a fake VCF file to give to R to do the Hierarchical Clustering then Sankoff algorithm on the Dendrogram

	:params possible_pairs: is a dictionary of position and the different amino acids at that position
	:params position_counts: dictionary of all positions the amino acid distribution at that count
	:params reads: dictionary of sample IDs and sequences
	'''

	header = ["Chromosome", "POS", "Indiv", "gt_PL", "gt_GT", "gt_GT_allele\n"]
	output_vcf = open("output_vcf.vcf" , "w")
	output_vcf.write("\t".join(header))
	# len_read = len(reads[list(reads.keys())[0]])
	for pos in possible_pairs.keys():
		for sample_id, read in reads.items():
				if read[pos] in possible_pairs[pos]:
					gt_GT = str(possible_pairs[pos].index(read[pos]))
				else:
					gt_GT = "NA"

				output_vcf.write("\t".join(["1", str(pos), str(sample_id), "0", gt_GT, ','.join(possible_pairs[pos] + ["\n"])]))


def generate_sequence_table(reads):
	'''
	outputs a TSV file for R to read with lines as samples and columns as positions in sequence 
	:params reads: dictionary of sample IDs and sequences
	'''
	output_table = open("input_protein_table.tsv", "w")
	# creating header for the table
	header = "\t".join([str(x) for x in range(len(reads[list(reads.keys())[1]]))])
	header = "samples_name\t" + header + "\n"
	output_table.write(header)
	for sample_id, read in reads.items():
		line = [str(sample_id)]
		line += [char for char in read]
		line[-1] = line[-1] + "\n"
		#pdb.set_trace()
		output_table.write("\t".join(line))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("You need to give the input reads file as first arguments and the percentage as second argument and 'dna' or 'aa' as third argument")
		sys.exit()
	if not os.path.exists(sys.argv[1]):
		print("The file {} does not exist".format(sys.argv[1]))
		sys.exit()

	reads = read_fasta(sys.argv[1])
	
	percentage = float(sys.argv[2])

	# generating the tsv table with each character in a column
	# needed for the R modules later
	generate_sequence_table(reads)

	if sys.argv[3] == "dna":
		position_counts = return_pos_counts(reads, percentage, dna=True)

	elif sys.argv[3] == "aa":
		position_counts = return_pos_counts(reads, percentage, dna=False)

	else:
		print("You need to give third argument as either 'dna' or 'aa'")
		sys.exit()
	
	total_count = len(reads)

	possible_pairs = dict()
	counter = 0
	for pos, aa_list in enumerate(position_counts):
		if len(aa_list) > 1:
			possible_pairs[pos] = [x.name for x in aa_list]  # this should be sorted according to counts
		# elif len(aa_list) == 1:
		# 	counter += 1
	# print("{} position had predominantly 1 protein from {} positions".format(counter, len(reads[list(reads.keys())[1]])))

	pairs_file = open("variant_pairs.txt", "w")
	pairs_file.write("first\tsecond\n")
	positions = list(possible_pairs.keys())
	# generating all possible pairs from position with variants
	# ignored position that are predominantly one amino acid
	for pair in itertools.combinations(positions, 2):
		pairs_file.write("{}\t{}\n".format(pair[0], pair[1]))
	pairs_file.close()

	generate_vcf_files(possible_pairs, position_counts, reads)
