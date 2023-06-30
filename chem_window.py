from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel, QWidget, QPushButton, QFileDialog, QMainWindow, QAction, QMenuBar, QInputDialog
from rdkit.Chem import Draw
from smiles_loader import SmilesLoader

class ChemWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.smiles_loader = SmilesLoader()
        self.current_mol_index = 0
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Chemical Structure Viewer')

        self.widget = QWidget()
        self.setCentralWidget(self.widget)

        self.layout = QVBoxLayout(self.widget)

        self.structure_layout = QHBoxLayout()
        self.label = QLabel()
        self.structure_layout.addWidget(self.label)
        self.layout.addLayout(self.structure_layout)

        self.button_layout = QHBoxLayout()
        self.prev_button = QPushButton('Previous', self)
        self.prev_button.clicked.connect(self.prev_structure)
        self.button_layout.addWidget(self.prev_button)

        self.next_button = QPushButton('Next', self)
        self.next_button.clicked.connect(self.next_structure)
        self.button_layout.addWidget(self.next_button)

        self.layout.addLayout(self.button_layout)

        self.load_csv_button = QPushButton('Load DB', self)
        self.load_csv_button.clicked.connect(self.load_local_csv)
        self.layout.addWidget(self.load_csv_button)

        self.create_menu()


    def create_menu(self):
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('File')

        loadAction = QAction('Load CSV', self)
        loadAction.triggered.connect(self.load_csv)
        fileMenu.addAction(loadAction)

        browseAction = QAction('Browse', self)
        browseAction.triggered.connect(self.browse_smiles)
        fileMenu.addAction(browseAction)
        

    def load_csv(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "CSV Files (*.csv);;All Files (*)", options=options)
        if fileName:
            self.smiles_loader.load_smiles(fileName)
            if len(self.smiles_loader) > 0:
                self.show_structure(0)

    def show_structure(self, index):
        if 0 <= index < len(self.smiles_loader):
            mol = self.smiles_loader.get_mol(index)
            pixmap = Draw.MolToQPixmap(mol)
            self.label.setPixmap(pixmap)
            self.label.repaint()  # force QLabel to update immediately
            self.current_mol_index = index
        else:
            print(f"Index out of range: {index}")

    def prev_structure(self):
        self.show_structure(self.current_mol_index - 1)

    def next_structure(self):
        self.show_structure(self.current_mol_index + 1)
        
    def browse_smiles(self):
        friendly_names_list = self.smiles_loader.get_all_friendly_names()
        friendly_name, okPressed = QInputDialog.getItem(self, "Browse","Chemical Structure:", friendly_names_list, 0, False)
        if okPressed and friendly_name:
            index = friendly_names_list.index(friendly_name)
            self.show_structure(index)
            
    def load_local_csv(self):
        # Attempt to load the local 'db.csv'
        self.smiles_loader.load_smiles('db.csv')
        if len(self.smiles_loader) > 0:
            self.show_structure(0)
