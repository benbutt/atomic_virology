import seq_tools

sars_cov_2 = seq_tools.fasta("./test_data/sequence.fasta")

sars_cov_2.find_orfs(200)