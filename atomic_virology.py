import seq_tools
import result_tools

# ## SEQ_TOOLS TESTING ###
# sars_cov_2 = seq_tools.fasta("./test_data/sequence.fasta")
# orfs = sars_cov_2.find_orfs(min_length=100)
# sars_cov_2.write_orf_fastas()

# ### RESULT_TOOLS TESTING ###
# test = result_tools.result("./test_data/ul51")
# test.get_results()
# test.get_plddts()
# test.plot_plddts()
# test.get_models()
# test.write_bfactors()
# test.get_msas()
# test.calculate_msa_depths()
# test.plot_msa_depth()

def main():
    alphafold_result = result_tools.result(args.result_dir)
    alphafold_result.get_results()
    alphafold_result.get_plddts()
    alphafold_result.plot_plddts()
    alphafold_result.get_models()
    alphafold_result.write_bfactors()
    alphafold_result.get_msas()
    alphafold_result.calculate_msa_depths()
    alphafold_result.plot_msa_depth()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Process raw AlphaFold2 output and build Jupyter notebook containing results")
    parser.add_argument("result_dir", type=str, help="Top-level directory containing AF2 output (usually target name)")
    args = parser.parse_args()

    main()