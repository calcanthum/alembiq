import os
import psycopg2
from dotenv import load_dotenv
from rdkit import Chem
from PyQt5.QtWidgets import QApplication, QMessageBox, QVBoxLayout, QWidget

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

def fetch_smiles(conn):
    cur = conn.cursor()
    cur.execute('SELECT id, smiles FROM molecules')
    smiles_data = cur.fetchall()
    cur.close()
    return smiles_data

def smiles_to_molblock(smiles_data):
    mol_blocks = []
    for data in smiles_data:
        try:
            mol = Chem.MolFromSmiles(data[1])
            mol_blocks.append((data[0], Chem.MolToMolBlock(mol)))  # save id and molblock as a tuple
        except Exception as e:
            return None, str(e)
    return mol_blocks, None

def update_molblocks(conn, mol_blocks):
    cur = conn.cursor()
    try:
        for block in mol_blocks:
            cur.execute('UPDATE molecules SET molblock = %s WHERE id = %s', (block[1], block[0]))
        conn.commit()
        return None
    except Exception as e:
        return str(e)
    finally:
        cur.close()

def create_message_box(error):
    app = QApplication([])
    msgBox = QMessageBox()
    if error:
        msgBox.setIcon(QMessageBox.Critical)
        msgBox.setText("Error")
        msgBox.setInformativeText(f"Failed to update molblock. Error: {error}")
    else:
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText("Success")
        msgBox.setInformativeText("Molblocks updated successfully.")
    msgBox.setStandardButtons(QMessageBox.Ok)
    msgBox.exec_()

def main():
    conn = connect_to_db()
    smiles_data = fetch_smiles(conn)
    mol_blocks, error = smiles_to_molblock(smiles_data)
    if error:
        create_message_box(error)
    else:
        error = update_molblocks(conn, mol_blocks)
        create_message_box(error)
    conn.close()

if __name__ == "__main__":
    main()
