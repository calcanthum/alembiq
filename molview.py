import os
from dotenv import load_dotenv
import psycopg2
from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtWidgets import (
    QApplication, QLabel, QVBoxLayout, QWidget, QPushButton, QFrame, QMenu, QFileDialog,
    QHBoxLayout, QTextEdit
)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
import db_core  # Import the db_core module


class MoleculeViewer(QWidget):
    def __init__(self):
        super().__init__()

        self.cur_molecule_index = 0  # Current index in the list of molecules
        self.molecules = []  # List of molecules

        self.layout = QVBoxLayout()  # GUI setup

        self.setup_middle_layout()
        self.setup_bottom_layout()

        self.setLayout(self.layout)

        self.load_molecules()
        self.display_molecule()

    def setup_middle_layout(self):
        self.middle_layout = QHBoxLayout()
        self.label = QLabel()
        self.label.setContextMenuPolicy(Qt.CustomContextMenu)
        self.label.customContextMenuRequested.connect(self.show_context_menu)
        self.middle_layout.addWidget(self.label)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.middle_layout.addWidget(self.info_text)

        self.layout.addLayout(self.middle_layout)

    def setup_bottom_layout(self):
        self.bottom_layout = QHBoxLayout()

        self.prev_button = QPushButton('Previous', self)
        self.prev_button.clicked.connect(self.prev_molecule)
        self.bottom_layout.addWidget(self.prev_button)

        self.next_button = QPushButton('Next', self)
        self.next_button.clicked.connect(self.next_molecule)
        self.bottom_layout.addWidget(self.next_button)

        self.layout.addLayout(self.bottom_layout)

    def load_molecules(self):
        connection = db_core.get_connection_from_pool()  # Get a connection from the pool
        if connection is None:
            return

        cursor = connection.cursor()

        try:
            cursor.execute("SELECT smiles, iupac FROM molecules")  # Select smiles and iupac from molecules
            rows = cursor.fetchall()
        except psycopg2.Error as e:
            print("Error fetching molecules: ", e)
            rows = []

        for row in rows:
            try:
                molecule = Chem.MolFromSmiles(row[0])
                if molecule is not None:
                    iupac_name = row[1] if row[1] is not None else "N/A"
                    molecule.SetProp("IUPACName", iupac_name)
                    self.molecules.append(molecule)
            except Exception as e:
                print(f"Error processing molecule: {e}")

        cursor.close()  # Close the cursor
        db_core.release_connection_to_pool(connection)  # Release the connection to the pool

    def display_molecule(self):
        # Display a molecule
        if self.molecules:
            molecule = self.molecules[self.cur_molecule_index]
            img = Draw.MolToQPixmap(molecule)
            self.label.setPixmap(img)
            self.display_molecule_info(molecule)
        else:
            self.label.setText("No molecules found in the database.")
            self.info_text.clear()

    def display_molecule_info(self, molecule):
        # Display the IUPAC name of the molecule
        iupac_name = molecule.GetProp("IUPACName")
        info = f"IUPAC Name: {iupac_name}"
        self.info_text.setText(info)

    def next_molecule(self):
        # Move to the next molecule
        if self.molecules and self.cur_molecule_index < len(self.molecules) - 1:
            self.cur_molecule_index += 1
        self.display_molecule()

    def prev_molecule(self):
        # Move to the previous molecule
        if self.molecules and self.cur_molecule_index > 0:
            self.cur_molecule_index -= 1
        self.display_molecule()

    def show_context_menu(self, pos):
        # Create a context menu
        menu = QMenu(self)

        # Add a "Save as PNG" option
        save_png_action = menu.addAction("Save as PNG")
        save_png_action.triggered.connect(self.save_image_as_png)

        # Display the context menu
        menu.exec_(self.label.mapToGlobal(pos))

    def save_image_as_png(self):
        # Save the current image as a PNG file
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