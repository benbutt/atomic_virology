import seq_tools

fasta = seq_tools.parse_fasta("test_data/sequence.fasta")
print(list(fasta))