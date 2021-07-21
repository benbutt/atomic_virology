from Bio import SeqIO

def parse_fasta(path):
    with open(path) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            yield record