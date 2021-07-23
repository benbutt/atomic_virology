from Bio import SeqIO, Seq
import re

class fasta:
    def __init__(self, path):
        self.path = path

        record = SeqIO.read(path, "fasta")
        self.id = record.id
        self.seq = record.seq
        self.length = len(self.seq)

    def find_orfs(self, min_length):
        """
        Indentify open reading frames in input fasta sequence
        Adapted from: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec384
        """
        print("Finding ORFs...")
        orfs = []

        for strand, nucs in [(+1, self.seq), (-1, self.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((self.length-frame) // 3)
                query_seq = str(nucs[frame:frame+length].translate())

                peptides = re.findall("M[^\*]+\*", query_seq)  
                peptides = [ Seq.Seq(peptide) for peptide in peptides if len(peptide) >= min_length ]
                orfs += peptides

        print(f"Found {len(orfs)} orfs!")
        self.orfs = orfs
