import os
import psycopg2
from dotenv import load_dotenv
from rdkit import Chem
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QTextEdit, QPushButton, QMessageBox

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

def fetch_smiles_and_molblocks(conn):
    cur = conn.cursor()
    cur.execute('SELECT id, smiles, molblock FROM molecules')
    data = cur.fetchall()
    cur.close()
    return data

def compare_smiles_and_molblocks(data):
    discrepancies = []
    for row in data:
        id_, smiles, molblock = row
        mol_from_molblock = Chem.MolFromMolBlock(molblock)
        smiles_from_molblock = Chem.MolToSmiles(mol_from_molblock)

        if smiles != smiles_from_molblock:
            discrepancies.append((id_, smiles, smiles_from_molblock))

    return discrepancies

def correct_discrepancies(conn, discrepancies):
    cur = conn.cursor()
    for id_, _, smiles_from_molblock in discrepancies:
        try:
            cur.execute('UPDATE molecules SET smiles = %s WHERE id = %s', (smiles_from_molblock, id_))
            conn.commit()
        except Exception as e:
            return str(e)
    cur.close()
    return None


def main():
    conn = connect_to_db()
    data = fetch_smiles_and_molblocks(conn)
    discrepancies = compare_smiles_and_molblocks(data)

    app = QApplication([])
    window = QWidget()
    layout = QVBoxLayout()
    text_edit = QTextEdit()
    button = QPushButton("Correct discrepancies")

    if discrepancies:
        details = ""
        for id_, smiles, smiles_from_molblock in discrepancies:
            details += f'ID: {id_}, Original SMILES: {smiles}, SMILES from molblock: {smiles_from_molblock}\n'
        text_edit.setText(details)
        button.clicked.connect(lambda: correct_discrepancies(conn, discrepancies))
    else:
        text_edit.setText("No discrepancies found.")
        button.setEnabled(False)

    layout.addWidget(text_edit)
    layout.addWidget(button)
    window.setLayout(layout)
    window.show()
    app.exec_()

    conn.close()

if __name__ == "__main__":
    main()
