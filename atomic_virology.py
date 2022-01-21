"""
Code for handling parsing AlphaFold 2 results, including structure alignment and plots
Ben Butt 2021
Updated: Jan 2022
"""

from result_tools import Result

def main():
    # Parse input dir
    alphafold_result = Result(args.result_dir)
    # Parse and tidy up results
    alphafold_result.get_results()
    # Extract scores from results pickle
    alphafold_result.get_scores()
    # Write pAE and pLDDT scores to .csv for convenience
    alphafold_result.write_scores()
    # Plot pLDDT scores
    alphafold_result.plot_plddt()
    # Plot pAE scores
    alphafold_result.plot_pae()
    # Get models and superimpose them on best model
    alphafold_result.get_models()
    alphafold_result.superimpose_models()
    # Parse MSAs and calculate per-position MSA depth
    alphafold_result.get_msas()
    alphafold_result.plot_msa_depth()
    # TODO: MSA handling for multimers

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Process raw AlphaFold2 output and build Jupyter notebook containing results")
    parser.add_argument("result_dir", type=str, help="Top-level directory containing AF2 output (usually target name)")
    args = parser.parse_args()

    main()