import seq_tools

test = seq_tools.fasta("./test_data/sequence.fasta")

test.find_orfs(10)

#print(test.orfs[0])