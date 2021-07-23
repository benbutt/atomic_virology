import seq_tools

sars_cov_2 = seq_tools.fasta("./test_data/sequence.fasta")
orfs = sars_cov_2.find_orfs(200)
sars_cov_2.write_orf_fastas("./test_data")