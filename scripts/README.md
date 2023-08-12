These scripts are used to further filter and analyze the results produced by the main tool in R.

 1. `calculate_pairwise_distances.py` takes a PDB file in and calculates all pairwise distances between all residues. This was used to get all the distances then look at the distances between the the pair of amino acids of interest
The output of this script is a TSV file with both residues id, name, and distance between them.

 2. `filter_p_value.py` takes a TSV file which is the output of the main tool's R script, this TSV file contain the pairs and their p-values. This script then filters and keeps the pairs with p-value up to a threshold given by the user as a second argument. Examle:
`python3 filter_p_value.py in_pairs.tsv 0.05` this will keep only pairs with a p-value up to 0.05

3. `multiple_testing_correction.py` this script corrects the p-value against multiple testing, multiple testing happens because a pair of positions can have several associated amnino acids and our tool will test these pair of positions several times for each combination of amino acids. It takes the TSV file output of the tool and outputs a similar one but with adjusted p-values.

4. `position_frequency.py` this script takes 3 arguments, an input MSA in FASTA format, first position as integer, and second position as integer, it then prints out the frequencies of each amino acid for those corrisponding positions (can be used just to check how many possible pairs exists between two positions and what is the most common amino acid for those positions).
An example from some MSA of a protein, after running the command `python3 position_frequence.py input_msa.fasta 10 20` to get the amino acids at positions 10 and 20.

```
First position has the following amino acids and their frequency:
W: 7112
L: 323
F: 4
K: 1
V: 2
G: 120
-: 13
R: 7
C: 1
Second position has the following amino acids and their frequency:
V: 7504
I: 71
F: 1
-: 3
G: 1
D: 1
C: 1
Y: 1

```

6. `two_pos_frequency.py` this script can be used with an MSA in FASTA format to count the specific combination of two positions. For example, in the previous example we looked at positions 10 and 20 and saw which amino acids and their counts, but if we wanted to know how many times does W at position 10 and V at position 20 come together, this script can do this. However, it does it by looking at 2 amino acids for each position and list all the possibilities.
Example with the same MSA of the previous example, running the following command `python3 two_pos_frequency.py input_msa.fasta 10 20 W L V I` will give the following output:

```
For W at 10 and V at 20 there are 7042
For W at 10 and I at 20 there are 64
For L at 10 and V at 20 there are 321
For L at 10 and I at 20 there are 1

```

7. `make_sif.py` this simple script takes the TSV file with pairs and thei p-values and produces a simple SIF file as a network to visualize for example with Cytoscape as a graph. This script take two inputs, a tsv file and a path to an output file to write to.