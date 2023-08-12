import os, sys


if len(sys.argv) < 3:
	print("You need to give two arguments <input_pairs.tsv> <output_graph.sif>")
	sys.exit()

if not os.path.exists(sys.argv[1]):
	print(f"The file {sys.argv[1]} does not exist")
	sys.exit()


in_file = sys.argv[1]
out_file = sys.argv[2]

pairs = set()
with open(in_file, "r") as infile:
	next(infile)  # skipping header
	for l in infile:
		l = l.strip().split()
		first = int(l[0])
		second = int(l[1])
		if first > second:
			pairs.add((first, second))
		else:
			pairs.add((second, first))

with open(out_file, "w") as outfile:
	for p in pairs:
		outfile.write(f"node{str(p[0])} pp node{str(p[1])}\n")
