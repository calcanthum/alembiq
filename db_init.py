import os
import sys
from dotenv import load_dotenv
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import AllChem
import psycopg2

# load environment variables
load_dotenv()

class DBInitWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Connect to the PostgreSQL database
        self.conn = psycopg2.connect(
            host=os.environ.get('DB_HOST'),
            database=os.environ.get('DB_NAME'),
            user=os.environ.get('DB_USER'),
            password=os.environ.get('DB_PASSWORD')
        )

        # Create a cursor object
        self.cursor = self.conn.cursor()

        # Delete the "molecules" table if it exists, then recreate it
        self.delete_table()
        self.create_table()

        # Set up the RDKit molecule
        smiles = 'CCO'
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates for the molecule

        # Insert the molecule into the database
        self.insert_molecule(mol)

        # Close the database connection
        self.cursor.close()
        self.conn.close()

        # Display a label
        label = QLabel('Database Initialized!', self)
        label.setAlignment(Qt.AlignCenter)
        self.setCentralWidget(label)

    def delete_table(self):
        # Prepare the SQL statement to drop the table
        statement = """
        DROP TABLE IF EXISTS molecules
        """

        # Execute the SQL statement
        self.cursor.execute(statement)
        self.conn.commit()

    def create_table(self):
        # Prepare the SQL statement for table creation
        statement = """
        CREATE TABLE IF NOT EXISTS molecules (
            id SERIAL PRIMARY KEY,
            smiles TEXT,
            molblock TEXT,
            inchi TEXT,
            iupac TEXT,
            friendly TEXT
        )
        """

        # Execute the SQL statement
        self.cursor.execute(statement)
        self.conn.commit()

    def insert_molecule(self, mol):
        # Generate the molecule's SMILES, molblock, InChI, IUPAC, and friendly name
        smiles = Chem.MolToSmiles(mol)
        molblock = Chem.MolToMolBlock(mol, includeStereo=True)
        inchi = Chem.MolToInchi(mol)
        iupac = "Ethanol"
        friendly = "Ethanol"

        # Remove the 3D coordinate tags from the molblock
        molblock_lines = molblock.split('\n')
        molblock_lines = [line for line in molblock_lines if not line.startswith('  3D')]

        # Join the molblock lines back together
        molblock = '\n'.join(molblock_lines)

        # Prepare the SQL statement for insertion
        statement = "INSERT INTO molecules (smiles, molblock, inchi, iupac, friendly) VALUES (%s, %s, %s, %s, %s)"
        values = (smiles, molblock, inchi, iupac, friendly)

        # Execute the SQL statement
        self.cursor.execute(statement, values)
        self.conn.commit()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = DBInitWindow()
    window.show()
    sys.exit(app.exec_())
