from PyQt5.QtWidgets import QLabel, QLineEdit, QPushButton, QVBoxLayout, QMessageBox, QWidget
from rdkit import Chem
from rdkit.Chem import Draw
import psycopg2
import os
from dotenv import load_dotenv

load_dotenv()

class SmilesDBEntry(QWidget):
    def __init__(self):
        super().__init__()

        # Connect to the PostgreSQL database
        self.conn = psycopg2.connect(
            dbname=os.getenv("DB_NAME"),
            user=os.getenv("DB_USER"),
            password=os.getenv("DB_PASSWORD"),
            host=os.getenv("DB_HOST"),
            port=os.getenv("DB_PORT")
        )

        self.initUI()

    def initUI(self):
        self.setWindowTitle('SMILES Database Entry')

        # QVBoxLayout for layout
        self.layout = QVBoxLayout()

        # Create a label, line edit (input field) and a submit button
        self.label = QLabel("Enter SMILES string:")
        self.smiles_input = QLineEdit()
        self.submit_button = QPushButton('Submit')

        # Connect the clicked signal of the button to the add_smiles method
        self.submit_button.clicked.connect(self.add_smiles)

        # Add widgets to layout
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.smiles_input)
        self.layout.addWidget(self.submit_button)

        # Set layout
        self.setLayout(self.layout)

    def add_smiles(self):
        smiles_string = self.smiles_input.text().strip()

        # Validate SMILES string using RDKit
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            QMessageBox.warning(self, 'Invalid SMILES', 'The entered SMILES string is not valid.')
            return

        # If SMILES string is valid, add it to the database
        cur = self.conn.cursor()
        cur.execute("INSERT INTO molecules (smiles) VALUES (%s)", (smiles_string,))
        self.conn.commit()

        # Show success message and clear input field
        QMessageBox.information(self, 'Success', 'SMILES string added to the database successfully.')
        self.smiles_input.clear()
