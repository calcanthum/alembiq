import sys
import os
from PyQt5.QtCore import QProcess
from PyQt5.QtWidgets import (
    QApplication, QVBoxLayout, QWidget, QPushButton, QMessageBox, QFrame
)
from db_init import DBInitWindow
from molview import MoleculeViewer
from db_add_molecule import MoleculeDBEntry
import db_smilesvsmolblock
from db_browser import DBBrowser
from moldraw import MoleculeDrawer

class Launcher(QWidget):
    def __init__(self):
        super().__init__()  # Call the constructor of the base class (QWidget)
        self.initUI()  # Initialize the User Interface
        self.windows = []  # List to hold references to open windows

    def initUI(self):
        self.setWindowTitle('Alembiq Launcher')  # Set the window title
        self.layout = QVBoxLayout()  # Create a vertical box layout
        
        # Buttons
        
        self.moldraw_button = QPushButton('Molecule Drawer', self) # Create a button for launching the molecule drawer
        self.moldraw_button.clicked.connect(self.moldraw) # Connect the button's click event to the moldraw method
        self.layout.addWidget(self.moldraw_button) # Add the button to the layout

        self.molview_button = QPushButton('Molecule Viewer', self) # Create a button for launching the molecule viewer
        self.molview_button.clicked.connect(self.molview) # Connect the button's click event to the molview method
        self.layout.addWidget(self.molview_button) # Add the button to the layout
        
        self.smiles_db_entry_button = QPushButton('Add molecule to DB', self) # Create a button for adding SMILES and InChI strings to the database
        self.smiles_db_entry_button.clicked.connect(self.smiles_db_entry) # Connect the button's click event to the smiles_db_entry method
        self.layout.addWidget(self.smiles_db_entry_button) # Add the button to the layout
                
        self.db_browser_button = QPushButton('DB Browser', self) 
        self.db_browser_button.clicked.connect(self.db_browser) 
        self.layout.addWidget(self.db_browser_button) 

        self.check_db_button = QPushButton('Check Database', self) # Create a button for checking the database
        self.check_db_button.clicked.connect(self.check_db) # Connect the button's click event to the check_db method
        self.layout.addWidget(self.check_db_button) # Add the button to the layout
                     
        self.line = QFrame() # Create a line to separate the sections
        self.line.setFrameShape(QFrame.HLine)
        self.line.setFrameShadow(QFrame.Sunken)
        self.layout.addWidget(self.line)
        
        self.init_db_button = QPushButton('\u26A0 Initialize Database', self) # Create a button for initializing the database
        self.init_db_button.clicked.connect(self.show_warning) # Connect the button's click event to the show_warning method
        self.layout.addWidget(self.init_db_button) # Add the button to the layout
        
        self.setLayout(self.layout) # Set the layout of the window to the layout we created

        
    def moldraw(self):
        self.molecule_drawer = MoleculeDrawer()  # Create a new instance of the MoleculeDrawer class
        self.molecule_drawer.show()  # Show the MoleculeDrawer
        self.windows.append(self.molecule_drawer)  # Add the MoleculeDrawer to our list of open windows

    def molview(self):
        self.molecule_viewer = MoleculeViewer()  # Create a new instance of the MoleculeViewer class
        self.molecule_viewer.show()  # Show the MoleculeViewer
        self.windows.append(self.molecule_viewer)  # Add the MoleculeViewer to our list of open windows

    def smiles_db_entry(self):
        self.smiles_db_entry = MoleculeDBEntry()  # Create a new instance of the MoleculeDBEntry class
        self.smiles_db_entry.show()  # Show the MoleculeDBEntry window
        self.windows.append(self.smiles_db_entry)  # Add the MoleculeDBEntry to our list of open windows
        
    def db_browser(self):
        self.db_browser = DBBrowser()  # Create a new instance of the MoleculeViewer class
        self.db_browser.show()  # Show the MoleculeViewer
        self.windows.append(self.db_browser)  # Add the MoleculeViewer to our list of open windows
        
    def check_db(self):
        # Launch the db_smilesvsmolblock.py script using QProcess
        self.process = QProcess(self)
        self.process.setWorkingDirectory(os.getcwd())
        self.process.errorOccurred.connect(self.process_error)  # Connect the errorOccurred signal to a custom slot
        self.process.start("python", ["db_smilesvsmolblock.py"])

    def process_error(self, error):
        print(f"Error occurred in process: {error}")

    def smiles_to_molblock(self):
        # Launch the db_smiles2molblock.py script using QProcess
        self.process = QProcess(self)
        self.process.setWorkingDirectory(os.getcwd())
        self.process.errorOccurred.connect(self.process_error)
        self.process.start("python", ["db_smiles2molblock.py"])
        print("Converting SMILES to Molblocks")

    def show_warning(self):
        choice = QMessageBox.warning(self, "Initialize Database",
                                     "This will wipe the 'molecules' table and the action is unrecoverable! Proceed?",
                                     QMessageBox.Yes | QMessageBox.No)
        if choice == QMessageBox.Yes:
            self.init_db()

    def init_db(self):
        self.db_window = DBInitWindow()  # Create a new instance of the DBInitWindow class
        self.db_window.show()  # Show the DBInitWindow
        self.windows.append(self.db_window)  # Add the DBInitWindow to our list of open windows

        
if __name__ == '__main__':
    app = QApplication(sys.argv)  # Create a new instance of the QApplication class
    ex = Launcher()  # Create a new instance of our Launcher class
    ex.show()  # Show the Launcher window
    sys.exit(app.exec_())  # Start the application's event loop
