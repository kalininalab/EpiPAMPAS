#!/usr/bin/python
import sys
import os
from Bio.PDB import *
from math import sqrt
import pdb


protein_letters_3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


if len(sys.argv) < 3:
    print("You need to give the following arguments <path_to_pdb.pdb> <path_to_output.tsv>")
    sys.exit()

if os.path.exists(sys.argv[1]):
    protein_name = sys.argv[1]
else:
    print("File {} does not exist".format(sys.argv[1]))
    sys.exit()

if os.path.exists(sys.argv[2]):
    print("The file {} will be overwritten".format(sys.argv[2]))

out_file = open(sys.argv[2], "w")


#Parse the pdb file
parser = PDBParser()
structure = parser.get_structure(protein_name, sys.argv[1])

residue_list = []

#Amino Acid Composition
for model in structure:
    # Iterate over the residues in the structure
    for residue in model.get_residues():
        # Check if residue isn't a HETATM
        if residue.get_id()[0] == " ":
            #residue_list.append(residue.get_resname())
            residue_list.append(residue)


distances = [["residue1 id", "residue1 name",
"resude2 id", "residue2 name", "distance"]]

for residue1 in residue_list:
    for residue2 in residue_list:
        if residue1 != residue2:

            distance = residue1['CA'] - residue2['CA']
            # pdb.set_trace()
            distances.append([str(residue1._id[1]), protein_letters_3to1[residue1.resname],
                 str(residue2._id[1]), protein_letters_3to1[residue2.resname],
                 str(distance)])

for l in distances:
    out_file.write("\t".join(l))
    out_file.write("\n")
