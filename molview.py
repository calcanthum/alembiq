import os
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtWidgets import QApplication, QLabel, QVBoxLayout, QWidget, QPushButton, QLineEdit, QFrame, QMenu, QFileDialog
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt
import psycopg2
import io

class MoleculeViewer(QWidget):
    def __init__(self):
        super().__init__()
        
        self.from_input = False # introduce a flag to keep track of whether the displayed molecule is from input or database
        self.cur_molecule_index = 0
        self.molecules = []

        self.layout = QVBoxLayout()
        self.label = QLabel()
        self.layout.addWidget(self.label)
        
        # Create a right-click menu
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)

        self.input_smiles = QLineEdit(self)
        self.input_smiles.setPlaceholderText("Enter SMILES")  # Add gray description text
        self.layout.addWidget(self.input_smiles)

        self.display_button = QPushButton('Display', self)
        self.display_button.clicked.connect(self.display_input_molecule)
        self.layout.addWidget(self.display_button)
        
        # Create a line to separate the sections
        self.line = QFrame() 
        self.line.setFrameShape(QFrame.HLine)
        self.line.setFrameShadow(QFrame.Sunken)
        self.layout.addWidget(self.line)

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
            self.next_button.setEnabled(False)  # Disable "Next" button
            self.prev_button.setEnabled(False)  # Disable "Previous" button
            return  # Exit the function early

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

    def display_input_molecule(self):
        smiles = self.input_smiles.text()
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            img = Draw.MolToQPixmap(molecule)
            self.label.setPixmap(img)
            self.from_input = True  # set the flag to True
        else:
            self.label.setText("Invalid SMILES string.")
        self.input_smiles.clear()

    def next_molecule(self):
        if self.from_input:  # if the current molecule is from input
            self.cur_molecule_index = 0  # reset the index to 0
            self.from_input = False  # reset the flag
        elif self.molecules and self.cur_molecule_index < len(self.molecules) - 1:
            self.cur_molecule_index += 1
        self.display_molecule()

    def prev_molecule(self):
        if self.from_input:  # if the current molecule is from input
            self.cur_molecule_index = len(self.molecules) - 1  # set the index to the last molecule
            self.from_input = False  # reset the flag
        elif self.molecules and self.cur_molecule_index > 0:
            self.cur_molecule_index -= 1
        self.display_molecule()

    def show_context_menu(self, pos):
        # Create the right-click menu
        menu = QMenu(self)

        # Add the "Save as PNG" option to the menu
        save_png_action = menu.addAction("Save as PNG")
        save_png_action.triggered.connect(self.save_image_as_png)

        # Show the menu at the clicked position
        menu.exec_(self.mapToGlobal(pos))

    def save_image_as_png(self):
        # Get the molecule image pixmap
        pixmap = self.label.pixmap()

        if pixmap is not None:
            # Get the file path using a file dialog or specify a default file path
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG (*.png)")

            if file_path:
                # Save the pixmap as a transparent PNG
                pixmap.save(file_path, "PNG", 0)  # 0 represents the quality (0-100)
        else:
            print("No image to save.")

if __name__ == "__main__":
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    app.exec_()