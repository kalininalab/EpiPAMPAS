import os
import sys


if __name__ == "__main__":
	aa = "ARNDBCEQZGHILKMFPSTWYV-"
	amino_acids = set(aa)
	for c in aa:
		amino_acids.add(c)
	
	if len(sys.argv) < 8:
		print("You need to give the following arguments: <input_msa.fasta> <first_pos> <second_pos> <first_pos_aa1> <first_pos_aa2> <second_pos_aa1> <second_pos_aa2>\n")
		print("For example if you want to investigate 2 positions in a protein (15, 67), and in position 15 the amino acids Z,A,R are usually present, and in position 67 H,I,K,L. You need to give 2 amino acids for each position and this script will return the frequency all the different combinations available and the number of sequences having them.\n")
		print("example input: input_table.tsv 15 67 Z R H K\n")
		print("The combination Z-H, Z-K, R-h, R-k and their ferquency will be returned")
		sys.exit()

	input_file = sys.argv[1]
	first_pos = int(sys.argv[2])
	second_pos = int(sys.argv[3])
	first_pos_aa = []
	first_pos_aa.append(sys.argv[4])
	first_pos_aa.append(sys.argv[5])
	second_pos_aa = []
	second_pos_aa.append(sys.argv[6])
	second_pos_aa.append(sys.argv[7])

	counters = [0,0,0,0]
	input_file = open(input_file, "r")
	for line in input_file.readlines():
		line = line.strip()
		if line[0] in amino_acids:
			if (line[first_pos] == first_pos_aa[0]) and (line[second_pos] == second_pos_aa[0]):
				counters[0] += 1
			elif (line[first_pos] == first_pos_aa[0]) and (line[second_pos] == second_pos_aa[1]):
				counters[1] += 1
			elif (line[first_pos] == first_pos_aa[1]) and (line[second_pos] == second_pos_aa[0]):
				counters[2] += 1
			elif (line[first_pos] == first_pos_aa[1]) and (line[second_pos] == second_pos_aa[1]):
				counters[3] += 1

	print("For {} at {} and {} at {} there are {}".format(first_pos_aa[0], first_pos, second_pos_aa[0], second_pos, counters[0]))
	print("For {} at {} and {} at {} there are {}".format(first_pos_aa[0], first_pos, second_pos_aa[1], second_pos, counters[1]))
	print("For {} at {} and {} at {} there are {}".format(first_pos_aa[1], first_pos, second_pos_aa[0], second_pos, counters[2]))
	print("For {} at {} and {} at {} there are {}".format(first_pos_aa[1], first_pos, second_pos_aa[1], second_pos, counters[3]))

	input_file.close()
