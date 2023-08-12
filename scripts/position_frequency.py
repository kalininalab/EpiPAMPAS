import os, gzip, sys


def read_fasta_gen(fasta_file_path):
    """
    A generator function that reads one read at a time
    Can be used for big FASTA files to not keep them in memory

    :param fasta_file_path: path to fasta file
    :yield: a tuple of sequence id and sequence
    """

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
                yield seq_name, "".join(seqs)
                seq_name = line[1:]
                seqs = []
            else:
                seq_name = line[1:]
        else:
            seqs.append(line)

    # last sequence
    if seqs:
        yield seq_name, "".join(seqs)


if len(sys.argv) < 4:
    print("You need to give the MSA in fasta format then first position then second position")
    sys.exit()

pos1 = int(sys.argv[2])
pos2 = int(sys.argv[3])

some_dict = {"pos1":{}, "pos2":{}}


counter = 0
for _, seq in read_fasta_gen(sys.argv[1]):
    if counter == 0:
        read_length = len(seq)
        counter += 1
    try:
        if pos1 < read_length:
            some_dict["pos1"][seq[pos1]] += 1
    except KeyError:
        if pos1 < read_length:
            some_dict["pos1"][seq[pos1]] = 1

    try:
        if pos2 < read_length:
            some_dict["pos2"][seq[pos2]] += 1
    except KeyError:
        if pos2 < read_length:
            some_dict["pos2"][seq[pos2]] = 1

print("First position has the following amino acids and their frequency:")
for k in some_dict['pos1'].keys():
    print(f"{k}: {some_dict['pos1'][k]}")

print("Second position has the following amino acids and their frequency:")
for k in some_dict['pos2'].keys():
    print(f"{k}: {some_dict['pos2'][k]}")