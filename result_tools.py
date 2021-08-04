import os
import shutil
from glob import glob
import pandas as pd
from matplotlib import pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

class result:
    """ 
    Class for handling raw AlphaFold outputs
    """
    def __init__(self, path):
        """
        Initialises result class with path to AF result directory
        """
        self.path = path

    def get_results(self):
        """
        Parses, stores and returns results pickle for each prediction
        """
        ## Move all raw output to a new subdirectory so it doesn't get in the way
        raw_files = glob(f"{self.path}/*") # Get all raw output files
        raw_dir = os.path.join(self.path, "raw_output") # Make a new subdirectory
        os.makedirs(raw_dir, exist_ok=True)     

        for raw_file in raw_files: # Iterate through raw files
            shutil.move(raw_file, raw_dir) # Move each file to new subdirectory

        ## Extract result pickle file for each prediction, parse and return
        result_paths = ( os.path.join(self.path, f"raw_output/result_model_{i}.pkl") for i in range(1,6) ) # Each path is generated on demand
        self.results = [ pd.read_pickle(result_path) for result_path in result_paths ] # Parse each results pickle using Pandas
        return self.results # Return full contents of results pickle as a list of dictionaries

    def get_plddts(self):
        """
        Extracts, stores and returns per-residue pLDDT scores from results pickle
        Requires: get_results()
        """
        ## Extract the pLDDT entry from each results dictionary
        self.plddts = [ result["plddt"] for result in self.results ]
        return self.plddts # Return pLDDT scores as a list of NumPy arrays

    def plot_plddts(self):
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
        # TODO: Add save as SVG for Inkscape
        plt.savefig(f"{plots_dir}/pLDDTs.png", dpi=600)
        # plt.show()

    def get_models(self):
        """
        Parses, stores and returns ranked models from results directory
        Requires: get_results()
        """
        parser = PDBParser() # Initialise a PDB parser object
        model_paths = ( os.path.join(self.path, f"raw_output/ranked_{i}.pdb") for i in range(5) ) # Generate the model paths on demand
        self.models = [ parser.get_structure(f"ranked_{i}", model_path) for i, model_path in enumerate(model_paths) ] # Parse each model and store as a list of models
        return self.models # Return a list of models

    def write_bfactors(self):
        """
        Writes out new PDB files with per-residue pLDDT scores in the B factor column
        Requires: get_results(), get_plddts(), get_models()
        """
        ## Set the b factor column for all atoms in each residue from the relevant pLDDT array
        for i, model in enumerate(self.models): # Iterate through the models
            for j, residue in enumerate(model.get_residues()): # Iterate through the residues from each model
                for atom in residue.get_atoms(): # Iterate through the atoms from each residue
                    atom.bfactor = self.plddts[i][j] # Set the B factor of each atom using the pLDDT array from get_plddts()

        ## Make a new subdirectory for storing model outputs
        models_dir = os.path.join(self.path, "models")
        os.makedirs(models_dir, exist_ok=True)
        
        ## Write out the new models with pLDDT score in B factor column
        io = PDBIO() # Initialise a PDBIO object for writing PDB files
        for i, model in enumerate(self.models): # Iterate through the models
            io.set_structure(model) # Get the current model
            io.save(f"{models_dir}/ranked_{i}_plddts.pdb") # Write each model in models subdirectory

    # TODO: Add MSA parsing including calculation and plotting of per-residue alignment depth