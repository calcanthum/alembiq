import sys
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QPushButton

class Launcher(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.windows = []

    def initUI(self):
        self.setWindowTitle('Program Launcher')

        self.layout = QVBoxLayout()
        
        self.structure_viewer_button = QPushButton('Structure Viewer', self)
        self.structure_viewer_button.clicked.connect(self.launch_structure_viewer)
        self.layout.addWidget(self.structure_viewer_button)
        
        self.setLayout(self.layout)

    def launch_structure_viewer(self):
        from chem_window import ChemWindow
        self.structure_viewer = ChemWindow()
        self.structure_viewer.show()
        self.windows.append(self.structure_viewer)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Launcher()
    ex.show()
    sys.exit(app.exec_())
