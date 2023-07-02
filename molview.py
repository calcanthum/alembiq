import os
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtWidgets import (
    QApplication, QLabel, QVBoxLayout, QWidget, QPushButton, QLineEdit, QFrame, QMenu, QFileDialog,
    QHBoxLayout, QTextEdit
)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
import psycopg2


class MoleculeViewer(QWidget):
    def __init__(self):
        super().__init__()

        self.from_input = False
        self.cur_molecule_index = 0
        self.molecules = []

        self.layout = QVBoxLayout()

        self.input_layout = QHBoxLayout()
        # Top 3rd layout
        self.layout.addLayout(self.input_layout)

        self.middle_layout = QHBoxLayout()

        self.label = QLabel()
        self.label.setContextMenuPolicy(Qt.CustomContextMenu)  # Set context menu policy
        self.label.customContextMenuRequested.connect(self.show_context_menu)  # Connect to the show_context_menu method
        self.middle_layout.addWidget(self.label)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.middle_layout.addWidget(self.info_text)

        self.layout.addLayout(self.middle_layout)

        self.line = QFrame()
        self.line.setFrameShape(QFrame.HLine)
        self.line.setFrameShadow(QFrame.Sunken)
        self.layout.addWidget(self.line)

        self.bottom_layout = QHBoxLayout()

        self.prev_button = QPushButton('Previous', self)
        self.prev_button.clicked.connect(self.prev_molecule)
        self.bottom_layout.addWidget(self.prev_button)

        self.next_button = QPushButton('Next', self)
        self.next_button.clicked.connect(self.next_molecule)
        self.bottom_layout.addWidget(self.next_button)

        self.layout.addLayout(self.bottom_layout)

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

        try:
            cursor = connection.cursor()
            cursor.execute("SELECT smiles, iupac FROM molecules")
            rows = cursor.fetchall()
        except psycopg2.DatabaseError as e:
            print(f"Database error: {e}")
            connection.close()
            self.next_button.setEnabled(False)
            self.prev_button.setEnabled(False)
            return

        for row in rows:
            try:
                molecule = Chem.MolFromSmiles(row[0])
                if molecule is not None:
                    iupac_name = row[1] if row[1] is not None else "N/A"
                    molecule.SetProp("IUPACName", iupac_name)  # Set the "iupac" value as a property
                    self.molecules.append(molecule)
            except Exception as e:
                print(f"Error processing molecule: {e}")

        connection.close()

    def display_molecule(self):
        if self.molecules:
            molecule = self.molecules[self.cur_molecule_index]
            img = Draw.MolToQPixmap(molecule)
            self.label.setPixmap(img)
            self.display_molecule_info(molecule)
        else:
            self.label.setText("No molecules found in the database.")
            self.info_text.clear()

    def display_molecule_info(self, molecule):
        iupac_name = molecule.GetProp("IUPACName")

        info = f"IUPAC Name: {iupac_name}"
        self.info_text.setText(info)

    def display_input_structure(self):
        input_str = self.input_field.text()
        try:
            if input_str.startswith('InChI='):
                molecule = Chem.MolFromInchi(input_str)
            else:
                molecule = Chem.MolFromSmiles(input_str)
        except Exception as e:
            print(f"Invalid chemical string: {e}")
            molecule = None

        if molecule is not None:
            try:
                img = Draw.MolToQPixmap(molecule)
                self.label.setPixmap(img)
                self.from_input = True
                self.display_molecule_info(molecule)
            except Exception as e:
                print(f"Error drawing molecule: {e}")
        else:
            self.label.setText("Invalid SMILES or InChI string.")
            self.info_text.clear()
        self.input_field.clear()

    def next_molecule(self):
        if self.from_input:
            self.cur_molecule_index = 0
            self.from_input = False
        elif self.molecules and self.cur_molecule_index < len(self.molecules) - 1:
            self.cur_molecule_index += 1
        self.display_molecule()

    def prev_molecule(self):
        if self.from_input:
            self.cur_molecule_index = len(self.molecules) - 1
            self.from_input = False
        elif self.molecules and self.cur_molecule_index > 0:
            self.cur_molecule_index -= 1
        self.display_molecule()

    def show_context_menu(self, pos):
        menu = QMenu(self)

        save_png_action = menu.addAction("Save as PNG")
        save_png_action.triggered.connect(self.save_image_as_png)

        menu.exec_(self.label.mapToGlobal(pos))  # Use label's mapToGlobal method

    def save_image_as_png(self):
        pixmap = self.label.pixmap()

        if pixmap is not None:
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG (*.png)")
            if file_path:
                try:
                    pixmap.save(file_path, "PNG", 0)
                except Exception as e:
                    print(f"Error saving image: {e}")
        else:
            print("No image to save.")

if __name__ == "__main__":
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    app.exec_()
