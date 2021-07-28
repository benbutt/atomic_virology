import pandas as pd


from matplotlib import pyplot as plt

class results:
    def __init__(self, path):
        self.path = path
        self.results = [ pd.read_pickle(f"{self.path}/result_model_{i}.pkl") for i in range(1,6) ]
        self.plddts = [ result["plddt"] for result in self.results ]

    def plot_plddts(self):
        fig, ax = plt.subplots()
        for i, plddts in enumerate(reversed(self.plddts)):
            plt.plot(range(1, len(plddts,)+1), plddts, label=f"Ranked_{4-i}.pdb")
        
        ax.set_xlabel("Residue number")
        ax.set_ylabel("pLDDT")
        plt.legend()
        plt.savefig("test_data/plots/pLDDTs.png", dpi=600)
        #plt.show()