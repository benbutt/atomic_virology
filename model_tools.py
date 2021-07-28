import Bio.PDB as pdb
import result_tools

class models:
    def __init__(self, path):
        self.path = path
        parser = pdb.PDBParser()
        self.models = [ parser.get_structure(f"ranked_{i}", f"{self.path}/ranked_{i}.pdb") for i in range(5) ]
    
    def write_plddts(self):
        for model in self.models:
            c_alphas = model.get_atoms().select("CA")
            for c_alpha in c_alphas:
                c_alpha.bfactor = 100
    # TODO : Read in pLDDT scores using result tools and stick in here