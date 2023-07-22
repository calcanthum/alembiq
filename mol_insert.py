from flask import jsonify, make_response
from db_core import get_connection_from_pool, release_connection_to_pool
from rdkit import Chem
import uuid

def insert_molecule(request):
    try:
        # Get the data from the POST request
        molecule_string = request.json.get('molecule_string', '')
        iupac = request.json.get('iupac', '')

        # Try to parse as InChI, if it fails, try SMILES
        mol = Chem.MolFromInchi(molecule_string)
        if mol is None:
            mol = Chem.MolFromSmiles(molecule_string)
        if mol is None:
            return make_response(jsonify({"error": "Invalid InChI or SMILES string."}), 400)

        # Generate the molecule's SMILES, molblock, InChI, and UUID
        smiles = Chem.MolToSmiles(mol)
        molblock = Chem.MolToMolBlock(mol, includeStereo=True)
        inchi = Chem.MolToInchi(mol)

        # Properly define and assign a value to the uuid variable
        generated_uuid = str(uuid.uuid4())

        conn = get_connection_from_pool()

        # Execute the insertion to the molecules table
        values = (generated_uuid, smiles, inchi, iupac)  # Use the generated_uuid
        statement = "INSERT INTO molecules (uuid, smiles, inchi, iupac) VALUES (%s, %s, %s, %s)"
        cursor = conn.cursor()
        cursor.execute(statement, values)
        conn.commit()

        # Get the molecule_id
        cursor.execute("SELECT id FROM molecules WHERE uuid = %s", (generated_uuid,))  # Use the generated_uuid
        molecule_id = cursor.fetchone()[0]

        # Insert molblock into the molblocks table
        values = (molecule_id, molblock)
        statement = "INSERT INTO molblocks (molecule_id, molblock) VALUES (%s, %s)"
        cursor.execute(statement, values)
        conn.commit()

        release_connection_to_pool(conn)

        return make_response(jsonify({"message": "Molecule inserted successfully."}), 200)
    except Exception as e:
        print(f"Error inserting molecule: {e}")
        return make_response(jsonify({"error": "Error inserting molecule."}), 500)
