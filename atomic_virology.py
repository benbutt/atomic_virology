import seq_tools
import result_tools

# ## SEQ_TOOLS TESTING ###
# sars_cov_2 = seq_tools.fasta("./test_data/sequence.fasta")
# orfs = sars_cov_2.find_orfs(min_length=100)
# sars_cov_2.write_orf_fastas()

### RESULT_TOOLS TESTING ###
test = result_tools.result("./test_data/ul51")
test.get_results()
test.get_plddts()
test.plot_plddts()
test.get_models()
test.write_bfactors()
test.get_msas()