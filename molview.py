import os
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtWidgets import QApplication, QLabel, QVBoxLayout, QWidget, QPushButton
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt
import psycopg2
import io

class MoleculeViewer(QWidget):
    def __init__(self):
        super().__init__()

        self.cur_molecule_index = 0
        self.molecules = []

        self.layout = QVBoxLayout()
        self.label = QLabel()
        self.layout.addWidget(self.label)

        self.next_button = QPushButton('Next', self)
        self.next_button.clicked.connect(self.next_molecule)
        self.layout.addWidget(self.next_button)

        self.prev_button = QPushButton('Previous', self)
        self.prev_button.clicked.connect(self.prev_molecule)
        self.layout.addWidget(self.prev_button)

        self.setLayout(self.layout)

        self.load_molecules()
        self.display_molecule()

    def load_molecules(self):
        load_dotenv()
        db_name = os.getenv('DB_NAME')
        db_user = os.getenv('DB_USER')
        db_password = os.getenv('DB_PASSWORD')
        db_host = os.getenv('DB_HOST')
        db_port = os.getenv('DB_PORT')

        connection = psycopg2.connect(
            dbname=db_name,
            user=db_user,
            password=db_password,
            host=db_host,
            port=db_port
        )

        cursor = connection.cursor()
        cursor.execute("SELECT smiles FROM molecules")
        rows = cursor.fetchall()

        for row in rows:
            molecule = Chem.MolFromSmiles(row[0])
            self.molecules.append(molecule)

        connection.close()

    def display_molecule(self):
        if self.molecules:
            molecule = self.molecules[self.cur_molecule_index]
            img = Draw.MolToQPixmap(molecule)
            self.label.setPixmap(img)
        else:
            self.label.setText("No molecules found in database.")


    def next_molecule(self):
        if self.molecules and self.cur_molecule_index < len(self.molecules) - 1:
            self.cur_molecule_index += 1
            self.display_molecule()

    def prev_molecule(self):
        if self.molecules and self.cur_molecule_index > 0:
            self.cur_molecule_index -= 1
            self.display_molecule()


if __name__ == "__main__":
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    app.exec_()