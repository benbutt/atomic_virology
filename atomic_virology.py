import seq_tools
import result_tools


"""
### SEQ_TOOLS TESTING ###
sars_cov_2 = seq_tools.fasta("./test_data/sequence.fasta")
orfs = sars_cov_2.find_orfs(100)
sars_cov_2.write_orf_fastas("./test_data/ORFs")
"""

"""
### RESULT_TOOLS TESTING ###
"""
test = result_tools.result("./test_data/ul51")
test.write_bfactors()

"""
### MODEL_TOOLS TESTING ###
ranked = model_tools.models("test_data/models)
"""