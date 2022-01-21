"""
Code for handling parsing AlphaFold 2 results, including structure alignment and plots
Ben Butt 2021
Updated: Jan 2022
"""

# Standard imports
import os
import subprocess
import shutil
import json
from glob import glob
from typing import List, Dict

# Third party imports
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

class result:
    """ 
    Class for handling raw AlphaFold outputs
    """
    def __init__(self, path: str) -> None:
        """
        Initialises result class with path to AF result directory
        """
        self.path = os.path.realpath(path)
        self.multimer = False
        print(f"Found result directory at {self.path}/")

    def get_results(self) -> List[Dict[str, np.ndarray]]:
        """
        Parses, stores and returns results pickle for each prediction
        """
        ## Move all raw output to a new subdirectory so it doesn't get in the way
        raw_files = glob(f"{self.path}/*") # Get all raw output files
        raw_dir = os.path.join(self.path, "raw_output") # Make a new subdirectory
        os.makedirs(raw_dir, exist_ok=True)     

        for raw_file in raw_files: # Iterate through raw files
            shutil.move(raw_file, raw_dir) # Move each file to new subdirectory
        print(f"Raw AlphaFold output moved to {raw_dir}/")

        ## Parse the rankings json to get the results in the correct order
        ranking_path = os.path.join(self.path, "raw_output", "ranking_debug.json")
        with open(ranking_path) as j:    
            ranking = json.load(j)["order"]

        ## Extract result pickle file for each prediction, parse and return
        result_paths = ( os.path.join(self.path, "raw_output", f"result_{rank}.pkl") for rank in ranking ) # Generate result paths on demand
        self.results = [ pd.read_pickle(result_path) for result_path in result_paths ] # Parse each results pickle using Pandas
        print(f"Parsed {len(self.results)} results files")
        return self.results # Return full contents of results pickle as a list of dictionaries        

    def get_scores(self) -> Dict[str, List[np.ndarray]]:
        """
        Extracts pLDDT, pAE and pTM/ipTM scores from results pickle
        Requires: get_results()

        Returns:
             Dict[str, List[np.ndarray]]: Dictionary containing lists of score arrays, keyed on score type
        """

        def _get_score(score_type: str) -> List[np.ndarray]:
            return [ result[score_type] for result in self.results ]

        if self.multimer:
            self.scores = {
                "pae": _get_score("predicted_aligned_error"), 
                "plddt": _get_score("plddt"), 
                "iptm": _get_score("iptm")
            }
        
        else:
            self.scores = {
                "pae": _get_score("predicted_aligned_error"), 
                "plddt": _get_score("plddt"), 
                "ptm": _get_score("ptm")
            }

        return self.scores

    def write_scores(self) -> None:
        """
        Writes pLDDT and pAE scores to .csv for further manipulation
        Requires: get_results(), get_scores
        """
        ## Make a new subdirectory for CSVs if it doesn't already exist
        csv_dir = os.path.join(self.path, "csv")
        os.makedirs(csv_dir, exist_ok=True)

        def _write_csv(score_type: str) -> None:
            scores = self.scores[score_type]
            
            for i in range(5):
                score = scores[i]
                csv_path = os.path.join(csv_dir, f"ranked_{i}_{score_type}.csv")
                pd.DataFrame(score).to_csv(csv_path)
        
        for score_type in ["pae", "plddt"]:
            _write_csv(score_type=score_type)


    def plot_plddt(self) -> Figure:
        """
        Plots per-residue pLDDT scores
        Requires: get_results(), get_scores()

        Returns: 
            plt.figure.Figure: Figure containing pLDDT plot
        """
        ## Make a new subdirectory for plots if it doesn't already exist
        plots_dir = os.path.join(self.path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        ## Initialise the plot object
        fig, ax = plt.subplots(constrained_layout=True)
        
        ## Set default font size to 8
        plt.rcParams.update({'font.size': 8})

        ## Do the plotting
        plddts = self.scores["plddt"] # Retrieve list of pLDDT arrays from scores dictionary
        ptms = self.scores["ptm"]
        res_nr = range(1, len(plddts[0])+1)
        for i, plddt in enumerate(plddts):
            ptm = self.scores
            ax.plot(res_nr, plddt, label=f"Model {i}, pTM: {ptms[i]:.2f}", zorder=0-i) # Plot each array against residue number, decreasing z-order for each plot (i.e. ranked_0 on top)
         
        ## Tidy up the plot
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        ax.set_ylim(0,100)
        ax.legend()
  
        # TODO: If multimer, mark chain termini on plot
        # if multimer:
        #     termini = [] # Calculate termini from sequence input?
        #     for terminus in termini:
        #         ax.axvline(terminus)

        ## Save the plot
        plt.savefig(f"{plots_dir}/pLDDTs.png", dpi=600)
        plt.savefig(f"{plots_dir}/pLDDTs.svg")
        print(f"Per-residue pLDDT plot saved to {plots_dir}/pLDDTs")

        return fig

    def plot_pae(self) -> Figure:
        """
        Plots per-residue pLDDT scores
        Requires: get_results(), get_plddts()

        Returns: 
            plt.figure.Figure: Figure containing pAE plot
        """
        ## Make a new subdirectory for plots if it doesn't already exist
        plots_dir = os.path.join(self.path, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        ## Initialise the plot object
        fig, axs = plt.subplots(nrows=1, ncols=5, sharey=True, constrained_layout=True)

        ## Set default font size to 8
        plt.rcParams.update({'font.size': 8})

        ## Do the plotting
        paes = self.scores["pae"]
        ptms = self.scores["ptm"]

        for i, pae in enumerate(paes):
            im = axs[i].imshow(pae, cmap="Greens_r", vmin=0, vmax=30) # assign im so we can add a cbar next to the last axis
            axs[i].set_title(f"Model {i}")

            ## Format the axes
            # Ticks
            axs[i].xaxis.set_major_locator(ticker.MultipleLocator(100))
            axs[i].xaxis.set_minor_locator(ticker.MultipleLocator(20))
            
            axs[i].yaxis.set_major_locator(ticker.MultipleLocator(100))
            axs[i].yaxis.set_minor_locator(ticker.MultipleLocator(20))
            
            # Labels
            axs[i].set_xlabel("pAE (Å)")
       
        # Format first axis
        axs[0].set_ylabel("Aligned residue")

        # Add colorbar
        plt.colorbar(im, shrink=0.2)

        ## Save the plot
        plt.savefig(f"{plots_dir}/pAE.png", dpi=75)
        plt.savefig(f"{plots_dir}/pAE.svg")
        print(f"Predicted aligned error plots saved to {plots_dir}/pAE")

        return fig

    def get_models(self) -> List[Structure]:
        """
        Parses, stores and returns ranked models from results directory
        Requires: get_results()
        """

        parser = PDBParser() # Initialise a PDB parser object
        model_paths = ( os.path.join(self.path, f"raw_output/ranked_{i}.pdb") for i in range(5) ) # Generate the model paths on demand
        self.models = [ parser.get_structure(f"ranked_{i}", model_path) for i, model_path in enumerate(model_paths) ] # Parse each model and store as a list of models
        print(f"Parsed {len(self.models)} models")
        return self.models # Return a list of models

    def superimpose_models(self) -> None:
        """
        Superimposes all models on best model (ranked_0)
        Requires: get_models()
        """
        ## Make a new subdirectory for aligned models if it doesn't already exist
        models_dir = os.path.join(self.path, "aligned_models")
        os.makedirs(models_dir, exist_ok=True)

        def _select_cas(model) -> List:
            return [ atom for atom in model.get_atoms() if atom.id == "CA" ]

        ## Initialise the classes we need to manipulate structures
        sup = Superimposer()
        io = PDBIO()

        ## Set the ref model and grab its CA atoms
        ref_model = self.models[0]
        test_models = self.models[1:]
        ref_CAs = _select_cas(model=ref_model)

        ## Superimpose each model on the best one
        for i, test_model in enumerate(test_models):
            test_CAs = _select_cas(test_model)
            sup.set_atoms(ref_CAs, test_CAs) # Calculate the rot/trans matrix
            print(f"Superimosed model {i+1}, RMSD={sup.rms:.2f}")
            sup.apply(test_model) # Apply the calculated matrix to the test model coords

            # Write out the aligned model coords
            model_path = os.path.join(models_dir, f"ranked_{i}_aligned.pdb")
            io.set_structure(test_model)
            io.save(model_path)
        
        ## Also write out the ref model since it is already "aligned"
        io.set_structure(ref_model)
        model_path = os.path.join(models_dir, "ranked_0_aligned.pdb")
        io.save(model_path)

    def convert_a3m(self, infile, outfile) -> None:
        """
        Converts alignment from .a3m to .sto format which is easier to parse
        Requires "reformat.pl" conversion script from hhsuite to convert a3m format
        """
        hhsuite_path = "/usr/local/xtal/hhsuite"
        reformat_script = os.path.join(hhsuite_path, "scripts/reformat.pl")

        subprocess.run([reformat_script, "a3m", "sto", infile, outfile],stdout=subprocess.DEVNULL)
        print(f"Converted {infile} to Stockholm format")

    def get_msas(self) -> Dict[str, MultipleSeqAlignment]:
        """
        Parses, stores and returns BFD, Mgnify, and Uniref90 MSAs from genetic searches
        Requires "reformat.pl" conversion script from hhsuite to convert a3m format
        Requires: get_results()
        """
        # TODO: If multimer, parse chain MSAs seperately 

        ## Locate MSA directory
        msa_dir = os.path.join(self.path, "raw_output/msas/")

        ## Convert bfd_uniclust_hits from .a3m to .sto format, which is easier to parse
        bfd_a3m_path = os.path.join(msa_dir, "bfd_uniclust_hits.a3m")
        bfd_sto_path = os.path.join(msa_dir, "bfd_uniclust_hits.sto")
        self.convert_a3m(infile=bfd_a3m_path, outfile=bfd_sto_path)

        ## Locate individual alignments
        bfd_path = os.path.join(msa_dir, "bfd_uniclust_hits.sto")
        mgnify_path = os.path.join(msa_dir, "mgnify_hits.sto")
        uniref90_path = os.path.join(msa_dir, "uniref90_hits.sto")

        ## Parse all of the alignments
        self.msas = {
            "bfd_hits" : AlignIO.read(bfd_path, "stockholm"),
            "mgnify_hits" : AlignIO.read(mgnify_path, "stockholm"),
            "uniref90_hits" : AlignIO.read(uniref90_path, "stockholm")
            }

        print(f"Parsed {len(self.msas)} MSAs")

        def _calculate_msa_depths(self) -> Dict[str, List[int]]:
            """
            Calculates, stores and returns number of non-zero entries of each column of each alignment from self.msas as a dictionary of lists
            Requires: get_results(), get_msas()
            """
            self.msa_depths = {} # Initialise an empty dictionary to store the msa depth
        
            ## Iterate over the alignments in self.msas, calculate and store per-column depth
            for alignment in self.msas:
                depths = [] # Initialise an empty list to keep track of the per-column depth
                alignment_array = np.array([list(record) for record in self.msas[alignment]]) # Convert each record in the alignment to a list, and store the list of lists as a numpy array

                for column in alignment_array.T != "-": # Transpose the array to iterate over columns, converting every non-gap entry to True (gaps = False)
                    depths.append(np.count_nonzero(column)) # Count the number of Trues (non-gap characters) in each row (column of orignial alignment)

                self.msa_depths[alignment] = depths # Store the list of per-position depths in the msa_depths dictionary initialised above
        
            print(f"Calculated MSA depth for {len(self.msas)} MSAs")

            return self.msa_depths # Return all per-column depth for each msa in self.msas as a dictionary of lists

        _calculate_msa_depths()

        return self.msas

    def plot_msa_depth(self) -> None:
        """
        Plots per-column depth of all alignments in self.msas
        Requires: get_results(), get_msas(), calculate_msa_depths()
        """
        ## Make a new subdirectory for plots if it doesn't already exist
        plots_dir = os.path.join(self.path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        ## Initialise the plot object
        fig, ax = plt.subplots()

        ## Do the plotting
        for msa_depth in self.msa_depths: # Iterate over the stored alignments
            depths = self.msa_depths[msa_depth] # Grab the per-residue depth
            ax.scatter(
                x = range(1, len(depths)+1),
                y = depths,
                s = 4,
                label=msa_depth
                )
            
        ## Tidy up the plot
        ax.set_xlabel("Alignment position")
        ax.set_ylabel("Depth (no. sequences)")
        plt.legend()
        plt.tight_layout()

        plt.savefig(f"{plots_dir}/msa_depth.png", dpi=600)
        plt.savefig(f"{plots_dir}/msa_depth.svg")
        print(f"MSA depth plot saved to {plots_dir}/msa_depth")
