import pandas as pd
from matplotlib import pyplot as plt
import Bio.PDB as pdb

class result:
    def __init__(self, path):
        self.path = path

    def get_results(self):
        self.results = [ pd.read_pickle(f"{self.path}/result_model_{i}.pkl") for i in range(1,6) ]
        return self.results

    def get_plddts(self):
        self.plddts = [ result["plddt"] for result in self.results ]
        return self.plddts

    def plot_plddts(self):
        self.get_plddts()
        fig, ax = plt.subplots()
        for i, plddts in enumerate(reversed(self.plddts)):
            plt.plot(range(1, len(plddts,)+1), plddts, label=f"Ranked_{4-i}.pdb")
        
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        plt.legend()
        plt.savefig("test_data/plots/pLDDTs.png", dpi=600)
        plt.show()

    def get_models(self):
        parser = pdb.PDBParser()
        self.models = [ parser.get_structure(f"ranked_{i}", f"{self.path}/ranked_{i}.pdb") for i in range(5) ]
        return self.models

    def write_bfactors(self):
        self.get_plddts()
        self.get_models()
        
        for i, model in enumerate(self.models):
            for j, residue in enumerate(model.get_residues()):
                for atom in residue.get_atoms():
                    atom.bfactor = self.plddts[i][j]