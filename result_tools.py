import os
import shutil
from glob import glob
import pandas as pd
from matplotlib import pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

class result:
    def __init__(self, path):
        self.path = path

    def get_results(self):
        raw_files = glob(f"{self.path}/*")

        raw_dir = os.path.join(self.path, "raw_output")
        os.makedirs(raw_dir, exist_ok=True)     

        for raw_file in raw_files:
            shutil.move(raw_file, raw_dir)

        result_paths = ( os.path.join(self.path, f"raw_output/result_model_{i}.pkl") for i in range(1,6) )
        self.results = [ pd.read_pickle(result_path) for result_path in result_paths ]
        return self.results

    def get_plddts(self):

        self.plddts = [ result["plddt"] for result in self.results ]
        return self.plddts

    def plot_plddts(self):

        plots_dir = os.path.join(self.path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        fig, ax = plt.subplots()

        for i, plddts in enumerate(self.plddts):
            plt.plot(range(1, len(plddts,)+1), plddts, label=f"Ranked_{i}.pdb", zorder=0-i)
        
        
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        plt.legend()

        plt.savefig(f"{plots_dir}/pLDDTs.png", dpi=600)
        # plt.show()

    def get_models(self):
        parser = PDBParser()
        model_paths = ( os.path.join(self.path, f"raw_output/ranked_{i}.pdb") for i in range(5) )
        self.models = [ parser.get_structure(f"ranked_{i}", model_path) for i, model_path in enumerate(model_paths) ]
        return self.models

    def write_bfactors(self):
            
        for i, model in enumerate(self.models):
            for j, residue in enumerate(model.get_residues()):
                for atom in residue.get_atoms():
                    atom.bfactor = self.plddts[i][j]

        models_dir = os.path.join(self.path, "models")
        os.makedirs(models_dir, exist_ok=True)
        
        io = PDBIO()
        for i, model in enumerate(self.models):
            io.set_structure(model)
            io.save(f"{models_dir}/ranked_{i}_plddts.pdb")