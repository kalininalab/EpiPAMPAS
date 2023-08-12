import os
import sys
import random
import pdb
import subprocess
import time


if len(sys.argv) < 4:
	print("You need to give the p-valus file (R output) first, then the pairwise distances "
		  "TSV file second and a p-value threshold to filter the first file third")
	sys.exit()

pvalues = sys.argv[1]
pairwise_distances = sys.argv[2]
threshold = float(sys.argv[3])


all_pairs = dict()
with open(pairwise_distances, "r") as in_file:
	next(in_file)  # skipping header

	for l in in_file:
		l = l.split("\t")
		pos1 = int(l[0])
		pos2 = int(l[2])
		dist = float(l[4])
		if pos1 > pos2:
			all_pairs[(pos1, pos2)] = dist
		else:
			all_pairs[(pos2, pos1)] = dist


all_epi_distances = []
with open(pvalues, "r") as in_file:
	next(in_file)  # skipping first line

	for l in in_file:
		l = l.split("\t")
		if (float(l[-1]) < threshold) or (float(l[-2]) < threshold):
			pos1 = int(l[0]) + 1
			pos2 = int(l[1]) + 1
			if pos1 > pos2:
				epistatic_pair = (pos1, pos2)

			else:
				epistatic_pair = (pos2, pos1)

			if epistatic_pair in all_pairs:
				all_epi_distances.append(all_pairs[epistatic_pair])

print("The number of pairs under the {} threshold are {}".format(threshold, len(all_epi_distances)))

# the wilcoxon test from scipy needs two lists the same size
# and R statistical packages are more trustworthy than the Python ones
# I tried sampling to make two lists the same size, but different samplings gave different results
# the R function always gives the same results
r_script = """
# Read in the data
options(warn=-1)

distances_file <- scan("tmp.txt", what="", sep="\n")

epi_distances = as.double(strsplit(distances_file[1], ", ")[[1]])
all_distances = as.double(strsplit(distances_file[2], ", ")[[1]])

w_test <- wilcox.test(epi_distances, all_distances)
print(paste0("The welcox test statistic is: ", w_test$statistic, " and the p-value is: ", w_test$p.value))

"""
out_r_script = open("tmp.r", "w")
out_r_script.write(r_script)
out_r_script.close()

out_file = open("tmp.txt", "w")
out_file.write(str(all_epi_distances)[1:-1] + "\n")
out_file.write(str(list(all_pairs.values()))[1:-1] + "\n")
out_file.close()
# print(wilcoxon(all_epi_distances_sampled, all_distances_sampled, zero_method="wilcox"))
subprocess.call("Rscript tmp.r", shell=True)
subprocess.call("rm tmp.r", shell=True)
subprocess.call("rm tmp.txt", shell=True)
