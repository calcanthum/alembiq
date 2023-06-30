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
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

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
            molblock TEXT
        )
        """

        # Execute the SQL statement
        self.cursor.execute(statement)
        self.conn.commit()

    def insert_molecule(self, mol):
        # Generate the molecule's SMILES and molblock
        smiles = Chem.MolToSmiles(mol)
        molblock = Chem.MolToMolBlock(mol)

        # Prepare the SQL statement for insertion
        statement = "INSERT INTO molecules (smiles, molblock) VALUES (%s, %s)"
        values = (smiles, molblock)

        # Execute the SQL statement
        self.cursor.execute(statement, values)
        self.conn.commit()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = DBInitWindow()
    window.show()
    sys.exit(app.exec_())
