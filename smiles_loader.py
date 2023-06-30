import pandas as pd
from rdkit import Chem

class SmilesLoader:
    def __init__(self):
        self.smiles_list = []
        self.friendly_names = []

    def load_smiles(self, filename):
        self.smiles_list = []
        self.friendly_names = []
        df = pd.read_csv(filename)
        self.smiles_list = df['smiles'].tolist()
        self.friendly_names = df['friendly'].tolist()

    def __len__(self):
        return len(self.smiles_list)

    def get_mol(self, index):
        if 0 <= index < len(self.smiles_list):
            smiles = self.smiles_list[index]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Failed to create a molecule from SMILES: {smiles}")
            return mol
        else:
            print(f"Index out of range: {index}")
            return None
        
    def get_all_smiles(self):
        return self.smiles_list.copy()
    
    def get_all_friendly_names(self):
        return self.friendly_names.copy()