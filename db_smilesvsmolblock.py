import os
import uuid
from dotenv import load_dotenv
from rdkit import Chem
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QTextEdit, QPushButton, QMessageBox
from functools import partial
import db_core
from db_core import get_connection_from_pool, release_connection_to_pool

db_core.init_connection_pool(minconn=1, maxconn=10)

load_dotenv()
print("Starting application...")
def fetch_smiles_molblocks_inchi_and_uuid(conn):
    cur = conn.cursor()
    cur.execute('SELECT id, smiles, molblock, inchi, uuid FROM molecules')
    data = cur.fetchall()
    cur.close()
    return data

def compare_smiles_molblocks_inchi_and_uuid(data):
    discrepancies = []
    errors = []
    missing_uuids = []
    for row in data:
        id_, smiles, molblock, inchi, molecule_uuid = row
        mol = None
        if molblock is not None:
            try:
                mol = Chem.MolFromMolBlock(molblock)
            except Exception as e:
                errors.append((id_, f"Malformed molblock: {str(e)}"))
        elif smiles is not None:
            try:
                mol = Chem.MolFromSmiles(smiles)
            except Exception as e:
                errors.append((id_, f"Malformed smiles: {str(e)}"))
        elif inchi is not None:
            try:
                mol = Chem.MolFromInchi(inchi)
            except Exception as e:
                errors.append((id_, f"Malformed inchi: {str(e)}"))

        if mol:
            smiles_from_mol = Chem.MolToSmiles(mol)
            inchi_from_mol = Chem.MolToInchi(mol)
            if not smiles or not inchi or not molblock:
                discrepancies.append((id_, smiles or smiles_from_mol, Chem.MolToMolBlock(mol), inchi or inchi_from_mol))
        else:
            errors.append((id_, "Missing molblock, smiles, and inchi. Cannot generate missing data."))

        if not molecule_uuid:
            missing_uuids.append(id_)

    return discrepancies, errors, missing_uuids

def correct_discrepancies(conn, discrepancies, missing_uuids):
    cur = conn.cursor()
    for id_, smiles_from_molblock, molblock_from_mol, inchi_from_molblock in discrepancies:
        try:
            cur.execute('UPDATE molecules SET smiles = %s, molblock = %s, inchi = %s WHERE id = %s', (smiles_from_molblock, molblock_from_mol, inchi_from_molblock, id_))
            conn.commit()
        except Exception as e:
            return f"Error: {str(e)}"
    
    for id_ in missing_uuids:
        try:
            new_uuid = str(uuid.uuid4())
            cur.execute('UPDATE molecules SET uuid = %s WHERE id = %s', (new_uuid, id_))
            conn.commit()
        except Exception as e:
            return f"Error: {str(e)}"
            
    cur.close()
    return "Successfully corrected discrepancies and added missing UUIDs."

def main():
    print("Inside main()")
    conn = get_connection_from_pool()  # fetch a connection from pool
    data = fetch_smiles_molblocks_inchi_and_uuid(conn)
    discrepancies, errors, missing_uuids = compare_smiles_molblocks_inchi_and_uuid(data)


    app = QApplication([])
    window = QWidget()
    layout = QVBoxLayout()
    text_edit = QTextEdit()
    button = QPushButton("Correct discrepancies and add missing UUIDs")

    def handle_button_click():
        result = correct_discrepancies(conn, discrepancies, missing_uuids)
        text_edit.setText(result)

    if discrepancies or errors or missing_uuids:
        details = ""
        for id_, smiles, molblock_from_mol, inchi in discrepancies:
            molblock_exists = "Present" if molblock_from_mol else "NULL"
            details += f'ID: {id_}, Original SMILES: {smiles if smiles else "NULL"}, Molblock from mol: {molblock_exists}, Original InChI: {inchi if inchi else "NULL"}\n'
        for id_ in missing_uuids:
            details += f'ID: {id_}, Missing UUID.\n'
        for id_, error in errors:
            details += f'ID: {id_}, Error: {error}\n'
        text_edit.setText(details)

        button.clicked.connect(handle_button_click)

    else:
        text_edit.setText("No discrepancies or errors found. No missing UUIDs.")
        button.setEnabled(False)

    layout.addWidget(text_edit)
    layout.addWidget(button)
    window.setLayout(layout)
    window.show()
    app.exec_()

    release_connection_to_pool(conn)  # release the connection back to pool


if __name__ == "__main__":
    main()
