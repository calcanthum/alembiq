from PyQt5.QtWidgets import QLabel, QLineEdit, QPushButton, QVBoxLayout, QMessageBox, QWidget
from rdkit import Chem
import psycopg2
import os
from dotenv import load_dotenv

load_dotenv()

class MoleculeDBEntry(QWidget):
    def __init__(self):
        super().__init__()

        # Connect to the PostgreSQL database
        load_dotenv()
        db_name = os.getenv('DB_NAME')
        db_user = os.getenv('DB_USER')
        db_password = os.getenv('DB_PASSWORD')
        db_host = os.getenv('DB_HOST')
        db_port = os.getenv('DB_PORT')

        if None in [db_name, db_user, db_password, db_host, db_port]:
            print("One or more required environment variables are missing.")
            self.next_button.setEnabled(False)
            self.prev_button.setEnabled(False)
            return

        try:
            connection = psycopg2.connect(
                dbname=db_name,
                user=db_user,
                password=db_password,
                host=db_host,
                port=db_port
            )
        except psycopg2.OperationalError as e:
            print(f"Could not connect to the database: {e}")
            self.next_button.setEnabled(False)
            self.prev_button.setEnabled(False)
            return

        self.initUI()

    def initUI(self):
        self.setWindowTitle('Molecule Database Entry')

        # QVBoxLayout for layout
        self.layout = QVBoxLayout()

        # Create labels, line edits (input fields), and a submit button
        self.label_smiles_inchi = QLabel("Enter SMILES or InChI string:")
        self.label_friendly = QLabel("Enter Friendly name (optional):")
        self.label_iupac = QLabel("Enter IUPAC name (optional):")
        self.molecule_input = QLineEdit()
        self.friendly_input = QLineEdit()
        self.iupac_input = QLineEdit()
        self.molecule_input.setPlaceholderText("Enter SMILES or InChI")
        self.friendly_input.setPlaceholderText("Enter Friendly name")
        self.iupac_input.setPlaceholderText("Enter IUPAC name")
        self.submit_button = QPushButton('Submit')

        # Connect the clicked signal of the button to the add_molecule method
        self.submit_button.clicked.connect(self.add_molecule)

        # Add widgets to layout
        self.layout.addWidget(self.label_smiles_inchi)
        self.layout.addWidget(self.molecule_input)
        self.layout.addWidget(self.label_friendly)
        self.layout.addWidget(self.friendly_input)
        self.layout.addWidget(self.label_iupac)
        self.layout.addWidget(self.iupac_input)
        self.layout.addWidget(self.submit_button)

        # Set layout
        self.setLayout(self.layout)

    def add_molecule(self):
        molecule_string = self.molecule_input.text().strip()
        friendly_name = self.friendly_input.text().strip()
        iupac_name = self.iupac_input.text().strip()

        # Validate SMILES or InChI string using RDKit
        mol_from_smiles = Chem.MolFromSmiles(molecule_string)
        mol_from_inchi = Chem.MolFromInchi(molecule_string)

        if mol_from_smiles is not None:
            molecule_type = "smiles"
            other_type = "inchi"
            other_string = Chem.MolToInchi(mol_from_smiles)
            molblock = Chem.MolToMolBlock(mol_from_smiles)
        elif mol_from_inchi is not None:
            molecule_type = "inchi"
            other_type = "smiles"
            other_string = Chem.MolToSmiles(mol_from_inchi)
            molblock = Chem.MolToMolBlock(mol_from_inchi)
        else:
            QMessageBox.warning(self, 'Invalid Input', 'The entered string is not a valid SMILES or InChI.')
            return

        # If input string is valid, add it to the database
        cur = self.conn.cursor()
        if friendly_name and iupac_name:
            cur.execute(f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, friendly, iupac) VALUES (%s, %s, %s, %s, %s)",
                        (molecule_string, other_string, molblock, friendly_name, iupac_name))
        elif friendly_name:
            cur.execute(f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, friendly) VALUES (%s, %s, %s, %s)",
                        (molecule_string, other_string, molblock, friendly_name))
        elif iupac_name:
            cur.execute(f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, iupac) VALUES (%s, %s, %s, %s)",
                        (molecule_string, other_string, molblock, iupac_name))
        else:
            cur.execute(f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock) VALUES (%s, %s, %s)",
                        (molecule_string, other_string, molblock))
        self.conn.commit()

        # Show success message and clear input fields
        QMessageBox.information(self, 'Success', f'{molecule_type.upper()} string added to the database successfully. The corresponding {other_type.upper()} and MolBlock string were also added.')
        self.molecule_input.clear()
        self.friendly_input.clear()
        self.iupac_input.clear()
