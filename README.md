
- [Introduction](#introduction)
  * [Input Schemes](#input-schemes)
    + [Variants and TF-TFBS pairs](#variants-and-tf-tfbs-pairs)
    + [Multiple sequence alignment (MSA)](#multiple-sequence-alignment-(MSA))

  * [Snakefile](#snakefile)

# Introduction
This tools allows to find epistasis in sequence sets based on a recurring mutation patterrns. For method's details, please see [LINK TO PAPER]. The initial version of the method was developed for epistatic events in transcription factors (TFs) and TF-binding sites (TFBSs).
The code can run on lists of sequence variants from VCF files, as well as on multiple protein or nucleotide sequence alignments and check for epistasis in all possible pairs of positions where variants exists.

## Input Schemes
There are 3 input schemes that the one can use to run this tool.

### Variants and TF-TFBS pairs
This input scheme uses the script [Epistasis with TF-TFBS](epistasis_with_TF_TFBS.R), which can be called using `Rscript` and takes as inputs:
- (`-v, --vcf`) a VCF file with variants generated from the different sequences of the organism studied.
- (`-p, --pairs`) a TSV with two columns for the TF-TFBS positions in the organism studied. If the tsv file have a header, then the argument `-l, --header` should be given as well
- (`-f, --flanking`) an integer representing the flanking regions before and after the TF and TFBS positions to be taken into consideration while looking for variants.
- (`-d, --dendrograms`) the output directory location
- (`-b, --binomial`) which binomial test should be used, wehre the default if this argument is not given is right-tailed binomial, user can specify other, e.g. left-tailed or two-tailed.
- (`-o, --out`) output file path and name

```shell script
Options:
	-v CHARACTER, --vcf=CHARACTER
		VCF file containing all genetic variants

	-p CHARACTER, --pairs=CHARACTER
		tab-delimited file containing the positions of a variant pair (one tab-delimited variant pair per line)

	-f FLANKING, --flanking=FLANKING
		The size of the flanking region to check starting from the TF variant site

	-l, --header
		does the variant pair file have a header line?

	-d CHARACTER, --dendrograms=CHARACTER
		directory in which the constructed dendrograms should be stored (dendrograms are not stored if this flag is not given)

	-b CHARACTER, --binomial=CHARACTER
		should the binomial test be two-tailed, left-tailed, or right-tailed (default: right-tailed)

	-o CHARACTER, --out=CHARACTER
		output file name

	-h, --help
		Show this help message and exit

```
### Multiple sequence alignment (MSA)
If you have seqeunces in an MSA in FASTA format, you first need to run the [preprocessing script](preprocessing_msa.py). This python script takes 3 arguments, the input FASTA file, a float as the percentage of samples that have a variant to be considered in the pairs generation (explained in the next paragraph) and the last argument is ("aa" or "dna") to specify what kind of inputs you are giving.
The preprocessing will generate three files
* A vcf files with the variants, their positions and which samples have those variants.
* An input sequence table as TSV, where the MSA is turned into a big tab-separated table
* A text files with possible pairs to examine.
The idea of the preprocessing is that every combination of pairs of positions in the MSA is checked, and at one position (a column in the MSA), if it's conserved up to the threshold given in the second argument (for example, if 0.05 was given, then if 95% of the sequences at that position have the same letter, that position will be considered as conserved) then this position will not be considered in checking for pairs. 

Once these output files are generated, the user can then run the [epistasis from generated files with preprocessing](epistasis_with_generated_preprocess.R). This R script takes the following arguments:
- (`-v, --vcf`) The vcf files generated with the preprocessing script.
- (`-p, --pairs`) The pairs table outputted by the preprocessing step.
- (`-t, --sequences_table`) The tab-delimited table outputted by the preprocessing script.
- (`-d, --dendrograms`) The output directory location
- (`-o, --out`) Output file path and name
- (`-s, --start`) From which pair to start checking in the pairs table given in option (`-p`) (explained in the next part)
- (`-l, --slice`) How many pairs to take after the start

```shell script

Options:
	-v CHARACTER, --vcf=CHARACTER
		VCF file containing all genetic variants

	-p PAIRS, --pairs=PAIRS
		tab-delimited file containing the positions of a variant pair (one tab-delimited variant pair per line)

	-t SEQUENCES_TABLE, --sequences_table=SEQUENCES_TABLE
		The protein sequences TSV generated with the python code

	-b CHARACTER, --binomial=CHARACTER
		should the binomial test be two-tailed, left-tailed, or right-tailed (default: right-tailed)

	-o OUTPUT, --output=OUTPUT
		You can specify the output table name

	-s START, --start=START
		Where to start calculating in the pairs table (in case there are too many pairs and you need run in parallel)

	-l SLICE, --slice=SLICE
		How many pairs to take after the start position)

	-h, --help
		Show this help message and exit

```

The idea of having a start and slice here, is that if there are too many possible pairs in the pairs table (e.g. 1000), you can then use the Snakefile to parellalize it by running 10 instances for example, and each slice being a hundred, so in the first one you run with start being 1, then 101, then 202 and so on, and the slice is 100. Simply the outputted tables can be concatinated as they have the same order and the same columns and column order.

<!-- ### Running the dendrogram algorithm on a VCF file with variants
Available options for [dendrogram based algorithm](dendrogram-based_algorithm.R)
```shell script
    Options:
	-v CHARACTER, --vcf=CHARACTER
		VCF file containing all genetic variants

	-p CHARACTER, --pairs=CHARACTER
		tab-delimited file containing the positions of a variant pair (one tab-delimited variant pair per line)

	-f FLANKING, --flanking=FLANKING
		The size of the flanking region to check starting from the TF variant site

	-l, --header
		does the variant pair file have a header line?

	-d CHARACTER, --dendrograms=CHARACTER
		directory in which the constructed dendrograms should be stored (dendrograms are not stored if this flag is not given)

	-b CHARACTER, --binomial=CHARACTER
		should the binomial test be two-tailed, left-tailed, or right-tailed (default: right-tailed)

	-o CHARACTER, --out=CHARACTER
		output file name

	-h, --help
		Show this help message and exit
```
### Running the dendrogram algorithm on aligned protein sequences
First, you need to use the preprocessing Python script on the aligned protein sequences to generate a "fake" VCF table file, a table for all possible pairs, and a matrix for all the sequences and position needed for filtering when running the R script.
The [proprocessing script](preprocess_proteins.py) only takes one argument which is the path to the aligned protein sequences file.
Available options for the [protein R script](protein_dendrogram_algorithm.R):
Options:
	-v CHARACTER, --vcf=CHARACTER
		VCF file containing all genetic variants

	-f FUNCTIONS, --functions=FUNCTIONS
		Path to the functions script

	-p PAIRS, --pairs=PAIRS
		tab-delimited file containing the positions of a variant pair (one tab-delimited variant pair per line)

	-r PROTEINS, --proteins=PROTEINS
		The protein sequences TSV generated with the python code

	-o OUTPUT, --output=OUTPUT
		You can specify the output table name

	-s START, --start=START
		Where to start calculating in the pairs table (in case there are too many pairs and you need run in parallel)

	-h, --help
		Show this help message and exit

For the first 3 options, the files are generated with the Python preprocessing script. The start is just an integer, because there are many combinations and it might take sometimes, to run in parallel you can give from which row in the available pairs table to start from and the script will run the next 1000 pairs and otuput the results in a table named according to the output option. For example, the variant pairs tables have 3320 lines, you can run 4 processes one with 1 as the start, another with 1001, another with 2002, and the last with 3003 which will run from 3003 to the end which 3320. This way one can run 4 processes in parallele then just merge the output tsv files using something like `sed 1d file.tsv  >> all.tsv` so you can remove the header and append the rest of the lines to `all.tsv`. -->

### Snakefile
I am using the Snakefile here to be able to run many combinations in parallel. The idea is that the number of variants pairs is very big as we are checking every possible pairwise combination in the aligned protiein sequences file. Therefore, we have thousands of combinations. The `-s START` option tells the Rscript where to start in the variant pairs table given in option `-p PAIRS` and `-l SLICE` tells Rscript how many pairs to take after the start point. Hence, the Snakefile can run many different start positions and each one would run from start to the next 1000 pairs.
