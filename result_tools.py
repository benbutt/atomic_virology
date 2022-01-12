import os
import subprocess
import shutil
from glob import glob
from typing import List, Dict

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

class result:
    """ 
    Class for handling raw AlphaFold outputs
    """
    def __init__(self, path: str) -> None:
        """
        Initialises result class with path to AF result directory
        """
        self.path = os.path.realpath(path)
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

        ## Extract result pickle file for each prediction, parse and return
        result_paths = sorted(glob(f"{self.path}/raw_output/result_model_*.pkl"))
        self.results = [ pd.read_pickle(result_path) for result_path in result_paths ] # Parse each results pickle using Pandas
        print(f"Parsed {len(self.results)} results files")
        return self.results # Return full contents of results pickle as a list of dictionaries

    def get_plddts(self) -> List[np.ndarray]:
        """
        Extracts, stores and returns per-residue pLDDT scores from results pickle
        Requires: get_results()
        """
        ## Extract the pLDDT entry from each results dictionary
        self.plddts = [ result["plddt"] for result in self.results ]
        print(f"Extracted pLDDT scores from {len(self.plddts)} results files")
        return self.plddts # Return pLDDT scores as a list of NumPy arrays

    def plot_plddts(self) -> None:
        """
        Plots per-residue pLDDT scores
        Requires: get_results(), get_plddts()
        """
        ## Make a new subdirectory for plots if it doesn't already exist
        plots_dir = os.path.join(self.path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        ## Initialise the plot object
        fig, ax = plt.subplots()

        ## Do the plotting
        for i, plddts in enumerate(self.plddts): # Iterate through the pLDDT arrays
            ax.plot(range(1, len(plddts,)+1), plddts, label=f"Ranked_{i}.pdb", zorder=0-i) # Plot each array against residue number, decreasing z-order for each plot (i.e. ranked_0 on top)

        ## Tidy up the plot
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        ax.set_ylim(0,100)
        plt.legend()
        plt.tight_layout()

        ## Save the plot
        plt.savefig(f"{plots_dir}/pLDDTs.png", dpi=600)
        plt.savefig(f"{plots_dir}/pLDDTs.svg")
        print(f"Per-residue pLDDT plot saved to {plots_dir}/pLDDTs")

        # TODO: If multimer, mark chain termini on plot

    def get_paes(self) -> List[np.ndarray]:
        """
        Extracts, stores and returns predicted aligned error scores from results pickle
        Requires: get_results()
        """
        ## Extract pAE scores from each results dictionary
        self.paes = [ result["predicted_aligned_error"] for result in self.results ]

        print(f"Extracted pAE scores from {len(self.paes)} results files")
        return self.paes # Return pAE scores as a list of NumPy arrays

    def plot_paes(self) -> None:
        """
        Plots per-residue pLDDT scores
        Requires: get_results(), get_plddts()
        """
        ## Make a new subdirectory for plots if it doesn't already exist
        plots_dir = os.path.join(self.path, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        ## Initialise the plot object
        fig, axs = plt.subplots(nrows=1, ncols=5, sharey=True, constrained_layout=True)

        for i, pae in enumerate(self.paes):
            axs[i].imshow(pae)
            axs[i].set_title(f"Model {i+1}")

            # TODO: Wrap up in function
            axs[i].xaxis.set_major_locator(ticker.MultipleLocator(50))
            axs[i].xaxis.set_minor_locator(ticker.MultipleLocator(10))
            
            axs[i].yaxis.set_major_locator(ticker.MultipleLocator(50))
            axs[i].yaxis.set_minor_locator(ticker.MultipleLocator(10))
            
        ## Save the plot
        plt.savefig(f"{plots_dir}/pAEs.png", dpi=600)
        plt.savefig(f"{plots_dir}/pAEs.svg")
        print(f"Predicted aligned error plots saved to {plots_dir}/pAEs")

    def get_models(self) -> List[Structure]:
        """
        Parses, stores and returns ranked models from results directory
        Requires: get_results()
        """

        #TODO: Superimpose models

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

        sup = Superimposer()
        io = PDBIO()

        ref_model = self.models[0]

        ref_CAs = [ atom for atom in ref_model.get_atoms() if atom.id == "CA" ]

        for i, test_model in enumerate(self.models[1:]):
            test_CAs = [ atom for atom in test_model.get_atoms() if atom.id == "CA" ]
            sup.set_atoms(ref_CAs, test_CAs)
            print(f"Superimosed model {i+2}, RMSD={sup.rms:.2f}")
            sup.apply(test_model)

            model_path = os.path.join(models_dir, f"ranked_{i+1}_aligned.pdb")
            io.set_structure(test_model)
            io.save(model_path)
            
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

        #TODO: If multimer, parse chain MSAs seperately

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

        return self.msas

    def calculate_msa_depths(self) -> Dict[str, List[int]]:
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

            # TODO: Sliding window averaging to smooth plot?
            
        ## Tidy up the plot
        ax.set_xlabel("Alignment position")
        ax.set_ylabel("Depth (no. sequences)")
        plt.legend()
        plt.tight_layout()

        plt.savefig(f"{plots_dir}/msa_depth.png", dpi=600)
        plt.savefig(f"{plots_dir}/msa_depth.svg")
        print(f"MSA depth plot saved to {plots_dir}/msa_depth")