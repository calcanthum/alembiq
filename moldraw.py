from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtWidgets import (
    QApplication, QLabel, QVBoxLayout, QWidget, QPushButton, QLineEdit, QHBoxLayout, QTextEdit
)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt


class MoleculeDrawer(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()

        self.input_layout = QHBoxLayout()
        self.input_field = QLineEdit(self)
        self.input_field.setPlaceholderText("Enter SMILES or InChI")
        self.input_layout.addWidget(self.input_field)

        self.display_button = QPushButton('Display', self)
        self.display_button.clicked.connect(self.display_input_structure)
        self.input_layout.addWidget(self.display_button)

        self.layout.addLayout(self.input_layout)

        self.middle_layout = QHBoxLayout()

        self.label = QLabel()
        self.middle_layout.addWidget(self.label)

        self.layout.addLayout(self.middle_layout)

        self.setLayout(self.layout)

    def display_molecule(self):
        if self.molecules:
            molecule = self.molecules[self.cur_molecule_index]
            img = Draw.MolToQPixmap(molecule)
            self.label.setPixmap(img)
            self.display_molecule_info(molecule)
        else:
            self.label.setText("No molecules found in the database.")
            self.info_text.clear()

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
    viewer = MoleculeDrawer()
    viewer.show()
    app.exec_()