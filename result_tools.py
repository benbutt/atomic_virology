import os
import shutil
from typing import List, Dict
from glob import glob
import pandas as pd
from matplotlib import pyplot as plt
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

class result:
    """ 
    Class for handling raw AlphaFold outputs
    """
    def __init__(self, path: str) -> None:
        """
        Initialises result class with path to AF result directory
        """
        self.path = path
        print(f"Found result directory at {path}/")

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
        result_paths = ( os.path.join(self.path, f"raw_output/result_model_{i}.pkl") for i in range(1,6) ) # Each path is generated on demand
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
        print(f"Extracted pLDDT scores from {len(self.results)} results files")
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
            plt.plot(range(1, len(plddts,)+1), plddts, label=f"Ranked_{i}.pdb", zorder=0-i) # Plot each array against residue number, decreasing z-order for each plot (i.e. ranked_0 on top)
        
        ## Tidy up the plot
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        plt.legend()
        plt.tight_layout()

        ## Save the plot
        plt.savefig(f"{plots_dir}/pLDDTs.png", dpi=600)
        plt.savefig(f"{plots_dir}/pLDDTs.svg")
        print(f"Per-residue pLDDT plot saved to {plots_dir}/pLDDTs")

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

    def write_bfactors(self) -> None:
        """
        Writes out new PDB files with per-residue pLDDT scores in the B factor column
        Requires: get_results(), get_plddts(), get_models()
        """
        ## Set the b factor column for all atoms in each residue from the relevant pLDDT array
        for i, model in enumerate(self.models): # Iterate through the models
            for j, residue in enumerate(model.get_residues()): # Iterate through the residues from each model
                for atom in residue.get_atoms(): # Iterate through the atoms from each residue
                    atom.bfactor = self.plddts[i][j] # Set the B factor of each atom using the pLDDT array from get_plddts()
        print(f"Wrote pLDDT scores to B factor column of {len(self.models)} models")
        ## Make a new subdirectory for storing model outputs
        models_dir = os.path.join(self.path, "models")
        os.makedirs(models_dir, exist_ok=True)
        
        ## Write out the new models with pLDDT score in B factor column
        io = PDBIO() # Initialise a PDBIO object for writing PDB files
        for i, model in enumerate(self.models): # Iterate through the models
            io.set_structure(model) # Get the current model
            io.save(f"{models_dir}/ranked_{i}_plddts.pdb") # Write each model in models subdirectory
        print(f"Saved {len(self.models)} models with pLDDT scores to {models_dir}")

    def get_msas(self) -> Dict[str, MultipleSeqAlignment]:
        """
        Parses, stores and returns Mgnify and Uniref90 MSAs from genetic searches
        Requires: get_results()
        """
        msa_dir = os.path.join(self.path, "raw_output/msas/")

        # TODO: Either parse proprietary a3m format or convert (.fasta, .sto?) first
        # bfd_path = os.path.join(msa_dir, "bfd_uniclust_hits.a3m")

        mgnify_path = os.path.join(msa_dir, "mgnify_hits.sto")
        uniref90_path = os.path.join(msa_dir, "uniref90_hits.sto")

        self.msas = {
            #"bfd_hits" : AlignIO.read(bfd_path, "a3m"), # Can't parse A3m currently
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
            ax.plot(range(1, len(depths)+1), depths, label=msa_depth) # Plot the per-residue depth against length of the alignment

        ## Tidy up the plot
        ax.set_xlabel("Alignment position")
        ax.set_ylabel("Depth (no. sequences)")
        plt.legend()
        plt.tight_layout()

        
        plt.savefig(f"{plots_dir}/msa_depth.png", dpi=600)
        plt.savefig(f"{plots_dir}/msa_depth.svg")
        print(f"MSA depth plot saved to {plots_dir}/msa_depth")