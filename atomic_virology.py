import seq_tools
import result_tools

"""
sars_cov_2 = seq_tools.fasta("./test_data/sequence.fasta")
orfs = sars_cov_2.find_orfs(100)
sars_cov_2.write_orf_fastas("./test_data")
"""

test = result_tools.results("test_data/results")
test.plot_plddts()