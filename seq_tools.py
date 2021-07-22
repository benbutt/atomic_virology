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
        for strand, nuc in [(+1, self.seq), (-1, self.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((self.length-frame) // 3)
                test_seq = str(nuc[frame:frame+length])

                print(re.findall("(?:(?!TAG|TAA|TGA)[A-SU-Z])+(?=TAG|TAA|TGA)", test_seq))



        for orf in orfs:
            print(orfs)
        

        print(f"Found {len(orfs)} ORFs in input fasta")
                    