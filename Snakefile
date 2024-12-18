# location of Rscript binary
Rscript_path = "/usr/bin/Rscript"

ranges = [1, 1001, 2002, 3003, 4004, 5005, 6006, 7007, 8008, 9009]

wildcard_constraints:
	range = "\d+"

rule all:
	input:
		expand("protein_pairs_output_{range}.tsv", range=ranges)

rule run_R:
	input:
		vcf = "some_path/output_vcf.vcf",  # output from preprocessing
		pairs = "some_path/variant_pairs.txt",  # output from preprocessing
		proteins = "some_path/input_protein_table.tsv",  # output from preprocessing

	output: "protein_pairs_output_{range}.tsv"

	shell:
		"/usr/bin/time -v {Rscript_path} protein_dendrogram_algorithm.R -v {input.vcf} -p {input.pairs} -r {input.proteins} -s {wildcards.range} -o {output} 2> {output}_stats.txt"
