from flask import jsonify, make_response
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdCoordGen

def get_molecule_structure(molecule_string):
    try:
        if molecule_string.startswith('InChI='):
            molecule = Chem.MolFromInchi(molecule_string)
        else:
            molecule = Chem.MolFromSmiles(molecule_string)
    except Exception as e:
        print(f"Invalid chemical string: {e}")
        return make_response(jsonify({"error": "Invalid chemical string."}), 400)

    try:
        rdCoordGen.AddCoords(molecule)
    except Exception as e:
        print(f"Error generating 3D coordinates: {e}")

    # Convert the molecule to PDB format
    AllChem.EmbedMolecule(molecule, useBasicKnowledge=True)
    pdb_data = Chem.MolToPDBBlock(molecule)

    # Get atomic symbols
    atom_data = [atom.GetSymbol() for atom in molecule.GetAtoms()]

    return jsonify({
        'pdb_data': pdb_data,
        'atom_data': atom_data
    })