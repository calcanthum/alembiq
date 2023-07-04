from PyQt5.QtWidgets import QLabel, QLineEdit, QPushButton, QVBoxLayout, QMessageBox, QWidget
from rdkit import Chem
import psycopg2
from dotenv import load_dotenv
import db_core  # Import the module itself, not the class
import uuid

load_dotenv()

class MoleculeDBEntry(QWidget):
    def __init__(self):
        super().__init__()

        try:
            db_core.init_connection_pool(minconn=1, maxconn=10)  # Call the function from the module
        except (ValueError, psycopg2.OperationalError) as e:
            print(e)
            return

        self.initUI()

    def calculate_uuid(self):
        return str(uuid.uuid4())

    def initUI(self):
        self.setWindowTitle('Molecule Database Entry')

        self.layout = QVBoxLayout()

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

        self.submit_button.clicked.connect(self.add_molecule)

        self.layout.addWidget(self.label_smiles_inchi)
        self.layout.addWidget(self.molecule_input)
        self.layout.addWidget(self.label_friendly)
        self.layout.addWidget(self.friendly_input)
        self.layout.addWidget(self.label_iupac)
        self.layout.addWidget(self.iupac_input)
        self.layout.addWidget(self.submit_button)

        self.setLayout(self.layout)

    def add_molecule(self):
        molecule_string = self.molecule_input.text().strip()
        friendly_name = self.friendly_input.text().strip()
        iupac_name = self.iupac_input.text().strip()

        if not molecule_string:
            QMessageBox.warning(self, 'Empty Input', 'Please enter a valid SMILES or InChI string.')
            return

        molecule_uuid = self.calculate_uuid()

        mol_from_smiles = Chem.MolFromSmiles(molecule_string)
        mol_from_inchi = Chem.MolFromInchi(molecule_string)

        conn = db_core.get_connection_from_pool()  # Getting connection from pool

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
            db_core.release_connection_to_pool(conn)  # Release connection if there's an error
            return

        result = db_core.add_molecule(conn, molecule_type, other_type, molecule_string, other_string, molblock,
                                      friendly_name, iupac_name, molecule_uuid)  # Using connection from pool

        db_core.release_connection_to_pool(conn)  # Release connection after use

        if not result:
            QMessageBox.warning(self, 'Error', 'An error occurred while adding the molecule to the database.')
            return

        QMessageBox.information(self, 'Success',
                                f'{molecule_type.upper()} string added to the database successfully. The corresponding {other_type.upper()} and MolBlock string were also added.')
        self.molecule_input.clear()
        self.friendly_input.clear()
        self.iupac_input.clear()
