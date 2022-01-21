import result_tools

def main():
    # Parse input dir
    alphafold_result = result_tools.result(args.result_dir)
    # Parse and tidy up results
    alphafold_result.get_results()
    # Extract and plot pAE scores
    alphafold_result.get_paes()
    alphafold_result.plot_paes()
    # Extract and plot pLDDT scores
    alphafold_result.get_plddts()
    alphafold_result.plot_plddts()
    # Get models and superimpose them on best model
    alphafold_result.get_models()
    alphafold_result.superimpose_models()

    # # TODO: MSA handling for multimers
    # # # alphafold_result.get_msas()
    # # # alphafold_result.calculate_msa_depths()
    # # # alphafold_result.plot_msa_depth()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Process raw AlphaFold2 output and build Jupyter notebook containing results")
    parser.add_argument("result_dir", type=str, help="Top-level directory containing AF2 output (usually target name)")
    args = parser.parse_args()

    main()