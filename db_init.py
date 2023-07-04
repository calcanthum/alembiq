import os
import sys
from dotenv import load_dotenv
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QMessageBox
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import AllChem
import psycopg2
from psycopg2 import sql
import uuid
import db_core  # Import the module itself

# load environment variables
load_dotenv()

class DBInitWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Get a database connection from the connection pool
        try:
            self.conn = db_core.get_connection_from_pool()
            # Create a cursor object
            self.cursor = self.conn.cursor()
        except (ValueError, psycopg2.OperationalError) as e:
            QMessageBox.critical(self, "Connection Failed", str(e))
            sys.exit(1)

        # Delete the "molecules" table if it exists, then recreate it
        self.delete_table()
        self.create_table()

        # Set up the RDKit molecule
        smiles = 'CCO'
        try:
            mol = Chem.MolFromSmiles(smiles)
            AllChem.EmbedMolecule(mol)  # Generate 3D coordinates for the molecule
        except:
            QMessageBox.critical(self, "RDKit Error", "Error processing the molecule.")
            sys.exit(1)

        # Insert the molecule into the database
        self.insert_molecule(mol)

        # Close the cursor and release the database connection back to the pool
        self.cursor.close()
        db_core.release_connection_to_pool(self.conn)

        # Display a label
        label = QLabel('Database Initialized!', self)
        label.setAlignment(Qt.AlignCenter)
        self.setCentralWidget(label)

    def delete_table(self):
        try:
            # Prepare the SQL statement to drop the table
            statement = """
            DROP TABLE IF EXISTS molecules
            """

            # Execute the SQL statement
            self.cursor.execute(statement)
            self.conn.commit()
        except psycopg2.Error as e:
            QMessageBox.critical(self, "Error Deleting Table", str(e))

    def create_table(self):
        try:
            # Prepare the SQL statement for table creation
            statement = """
            CREATE TABLE IF NOT EXISTS molecules (
                id SERIAL PRIMARY KEY,
                uuid UUID,
                molblock TEXT,
                smiles TEXT,
                inchi TEXT,
                iupac TEXT,
                friendly TEXT
            )
            """

            # Execute the SQL statement
            self.cursor.execute(statement)
            self.conn.commit()
        except psycopg2.Error as e:
            QMessageBox.critical(self, "Error Creating Table", str(e))

    def calculate_uuid(self):
        # Calculate UUID using uuid.uuid4()
        return str(uuid.uuid4())

    def insert_molecule(self, mol):
        try:
            # Generate the molecule's SMILES, molblock, InChI, IUPAC, and friendly name
            smiles = Chem.MolToSmiles(mol)
            molblock = Chem.MolToMolBlock(mol, includeStereo=True)
            inchi = Chem.MolToInchi(mol)
            iupac = "Ethanol"
            friendly = "Ethanol"

            # Prepare the SQL statement for insertion
            statement = sql.SQL("INSERT INTO molecules (smiles, molblock, inchi, iupac, friendly, uuid) VALUES (%s, %s, %s, %s, %s, %s)")
            values = (smiles, molblock, inchi, iupac, friendly, self.calculate_uuid())

            # Execute the SQL statement
            self.cursor.execute(statement, values)
            self.conn.commit()
        except psycopg2.Error as e:
            QMessageBox.critical(self, "Error Inserting Molecule", str(e))

if __name__ == '__main__':
    try:
        db_core.init_connection_pool(minconn=1, maxconn=10)  # Initiate the connection pool
        app = QApplication(sys.argv)
        window = DBInitWindow()
        window.show()
        sys.exit(app.exec_())
    finally:
        db_core.close_all_conns()  # Close all connections when the app is shutting down
