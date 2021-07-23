from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
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
                peptides = [ Seq(peptide) for peptide in peptides if len(peptide) >= min_length ]
                orfs += peptides

        print(f"Found {len(orfs)} orfs!")
        self.orfs = orfs
        return self.orfs

    def write_orf_fastas(self, path):
        for i, orf in enumerate(self.orfs):
            record = SeqRecord(
                orf,
                id=f"ORF_{i}",
                name=f"ORF_{i}")
            SeqIO.write(record, f"{path}/ORF_{i}.fasta", "fasta")
        print(f"Wrote {len(self.orfs)} ORFs to .fasta")