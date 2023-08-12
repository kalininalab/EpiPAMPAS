import os
import sys


if len(sys.argv) < 3:
	print("You need to give two arguments <in_tsv_file> and <cutoff_pvalue_threshold>")
	print("The TSV file is the same as the output of the main R script")
	print("Example: python3 filter_p_values.py test.tsv 0.05")
	sys.exit()


if os.path.exists(sys.argv[1]):
	in_file = sys.argv[1]
else:
	print(f"The path {sys.argv[1]} does not exist")
	sys.exit()

threshold = float(sys.argv[2])
# only tsv files
# files = [x for x in [f for f in os.listdir(in_directory) if os.path.isfile(f)] if x.endswith("t.tsv")]


threshold = float(sys.argv[1])

# for f in files:
out_file = open(in_file.split(os.sep)[-1].split(".")[1] + "_filtered.tsv", "w")
with open(in_file, "r") as infile:

	for l in infile:
		try:
			if (float(l.split("\t")[-1]) < threshold) or (float(l.split("\t")[-2]) <= threshold):
				out_file.write(l)
		except:
			# for the header line
			out_file.write(l)

out_file.close()
