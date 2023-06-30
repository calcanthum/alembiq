import sys
from PyQt5.QtCore import QProcess
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QPushButton
from chem_window import ChemWindow
from db_init import DBInitWindow
from molview import MoleculeViewer
from db_add_smiles import SmilesDBEntry

class Launcher(QWidget):
    def __init__(self):
        super().__init__()  # Call the constructor of the base class (QWidget)
        self.initUI()  # Initialize the User Interface
        self.windows = []  # List to hold references to open windows

    def initUI(self):
        self.setWindowTitle('Program Launcher')  # Set the window title
        self.layout = QVBoxLayout()  # Create a vertical box layout
        
        # Create a button for launching the molecule viewer
        self.molview_button = QPushButton('Molecule Viewer', self)
        # Connect the button's click event to the molview method
        self.molview_button.clicked.connect(self.molview)
        # Add the button to the layout
        self.layout.addWidget(self.molview_button)
        
        # Create a button for adding SMILES strings to the database
        self.smiles_db_entry_button = QPushButton('Add SMILES to DB', self)
        # Connect the button's click event to the smiles_db_entry method
        self.smiles_db_entry_button.clicked.connect(self.smiles_db_entry)
        # Add the button to the layout
        self.layout.addWidget(self.smiles_db_entry_button)
        
        # Create a button for launching the structure viewer
        self.structure_viewer_button = QPushButton('Structure Viewer', self)
        # Connect the button's click event to the launch_structure_viewer method
        self.structure_viewer_button.clicked.connect(self.launch_structure_viewer)
        # Add the button to the layout
        self.layout.addWidget(self.structure_viewer_button)
        
        # Create a button for initializing the database
        self.init_db_button = QPushButton('Initialize Database', self)
        # Connect the button's click event to the init_db method
        self.init_db_button.clicked.connect(self.init_db)
        # Add the button to the layout
        self.layout.addWidget(self.init_db_button)
        
        # Set the layout of the window to the layout we created
        self.setLayout(self.layout)

    def molview(self):
        self.molecule_viewer = MoleculeViewer()  # Create a new instance of the MoleculeViewer class
        self.molecule_viewer.show()  # Show the MoleculeViewer
        self.windows.append(self.molecule_viewer)  # Add the MoleculeViewer to our list of open windows

    def smiles_db_entry(self):
        self.smiles_db_entry = SmilesDBEntry()  # Create a new instance of the SmilesDBEntry class
        self.smiles_db_entry.show()  # Show the SmilesDBEntry window
        self.windows.append(self.smiles_db_entry)  # Add the SmilesDBEntry to our list of open windows

        
    def launch_structure_viewer(self):
        self.structure_viewer = ChemWindow()  # Create a new instance of the ChemWindow class
        self.structure_viewer.show()  # Show the ChemWindow
        self.windows.append(self.structure_viewer)  # Add the ChemWindow to our list of open windows

    def init_db(self):
        self.db_window = DBInitWindow()  # Create a new instance of the DBInitWindow class
        self.db_window.show()  # Show the DBInitWindow
        self.windows.append(self.db_window)  # Add the DBInitWindow to our list of open windows


if __name__ == '__main__':
    app = QApplication(sys.argv)  # Create a new instance of the QApplication class
    ex = Launcher()  # Create a new instance of our Launcher class
    ex.show()  # Show the Launcher window
    sys.exit(app.exec_())  # Start the application's event loop
