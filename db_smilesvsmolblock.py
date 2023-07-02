import os
import psycopg2
from dotenv import load_dotenv
from rdkit import Chem
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QTextEdit, QPushButton, QMessageBox
from functools import partial

load_dotenv()

DB_NAME = os.getenv("DB_NAME")
DB_USER = os.getenv("DB_USER")
DB_PASSWORD = os.getenv("DB_PASSWORD")
DB_HOST = os.getenv("DB_HOST")
DB_PORT = os.getenv("DB_PORT")

def connect_to_db():
    conn = psycopg2.connect(
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD,
        host=DB_HOST,
        port=DB_PORT
    )
    return conn

def fetch_smiles_molblocks_and_inchi(conn):
    cur = conn.cursor()
    cur.execute('SELECT id, smiles, molblock, inchi FROM molecules')
    data = cur.fetchall()
    cur.close()
    return data

def compare_smiles_molblocks_and_inchi(data):
    discrepancies = []
    errors = []
    for row in data:
        id_, smiles, molblock, inchi = row
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
    return discrepancies, errors


def correct_discrepancies(conn, discrepancies):
    cur = conn.cursor()
    for id_, smiles_from_molblock, molblock_from_mol, inchi_from_molblock in discrepancies:
        try:
            cur.execute('UPDATE molecules SET smiles = %s, molblock = %s, inchi = %s WHERE id = %s', (smiles_from_molblock, molblock_from_mol, inchi_from_molblock, id_))
            conn.commit()
        except Exception as e:
            return f"Error: {str(e)}"
    cur.close()
    return "Successfully corrected discrepancies."


def main():
    conn = connect_to_db()
    data = fetch_smiles_molblocks_and_inchi(conn)
    discrepancies, errors = compare_smiles_molblocks_and_inchi(data)

    app = QApplication([])
    window = QWidget()
    layout = QVBoxLayout()
    text_edit = QTextEdit()
    button = QPushButton("Correct discrepancies")

    def handle_button_click():
        result = correct_discrepancies(conn, discrepancies)
        text_edit.setText(result)

    if discrepancies or errors:
        details = ""
        for id_, smiles, molblock_from_mol, inchi in discrepancies:
            molblock_exists = "Present" if molblock_from_mol else "NULL"
            details += f'ID: {id_}, Original SMILES: {smiles if smiles else "NULL"}, Molblock from mol: {molblock_exists}, Original InChI: {inchi if inchi else "NULL"}\n'
        for id_, error in errors:
            details += f'ID: {id_}, Error: {error}\n'
        text_edit.setText(details)

        # Change this line to use the new handle_button_click function
        button.clicked.connect(handle_button_click)

    else:
        text_edit.setText("No discrepancies or errors found.")
        button.setEnabled(False)

    layout.addWidget(text_edit)
    layout.addWidget(button)
    window.setLayout(layout)
    window.show()
    app.exec_()

    conn.close()

if __name__ == "__main__":
    main()
