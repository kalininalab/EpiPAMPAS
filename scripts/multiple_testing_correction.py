import sys
import os
import pdb


if len(sys.argv) < 2:
	print("You need to give an input TSV file")
	sys.exit()

if os.path.exists(sys.argv[1]):
	in_file = sys.argv[1]
else:
	print(f"The file {sys.argv[1]} does not exist")
	sys.exit()


counts = dict()

out_file = open(in_file.split(os.sep)[-1].split(".")[0] + "_corrected_multiple_testing.tsv", "w")
# looping through the files and recording the positions that were multiply tested to be corrected
with open(in_file, "r") as infile:
	
	for l in infile:
		l = l.strip().split("\t")
		# header
		if l[0] == "pos1":
			continue
		else:
			# storing the pair and counting
			if int(l[0]) > int(l[1]):
				pair = (int(l[0]), int(l[1]))
			else:
				pair = (int(l[1]), int(l[0]))

			if pair in counts:
				counts[pair] += 1
			else:
				counts[pair] = 1


# correcting the p-values by multiplying by the number of test
with open(in_file, "r") as infile:
	
	for l in infile:
		if l.startswith("p"):
			out_file.write(l)
		else:

			l = l.strip().split("\t")
			# storing the pair and counting
			if int(l[0]) > int(l[1]):
				pair = (int(l[0]), int(l[1]))
			else:
				pair = (int(l[1]), int(l[0]))

			# correcting p-value by multiplying with number of test
			# then capping the p-value at 1
			l[-2] = float(l[-2]) * counts[pair]
			if l[-2] > 1:
				l[-2] = 1
			l[-1] = float(l[-1]) * counts[pair]
			if l[-1] > 1:
				l[-1] = 1

			out_file.write("\t".join([str(x) for x in l]))
			out_file.write("\n")

out_file.close()
