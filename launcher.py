import sys
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QPushButton
from chem_window import ChemWindow

class Launcher(QWidget):
    def __init__(self):
        super().__init__()  # Call the constructor of the base class (QWidget)
        self.initUI()  # Initialize the User Interface
        self.windows = []  # List to hold references to open windows

    def initUI(self):
        self.setWindowTitle('Program Launcher')  # Set the window title
        self.layout = QVBoxLayout()  # Create a vertical box layout
        
        # Create a button for launching the structure viewer
        self.structure_viewer_button = QPushButton('Structure Viewer', self)
        
        # Connect the button's click event to the launch_structure_viewer method
        self.structure_viewer_button.clicked.connect(self.launch_structure_viewer)
        
        # Add the button to the layout
        self.layout.addWidget(self.structure_viewer_button)
        
        # Set the layout of the window to the layout we created
        self.setLayout(self.layout)

    def launch_structure_viewer(self):
        self.structure_viewer = ChemWindow()  # Create a new instance of the ChemWindow class
        self.structure_viewer.show()  # Show the ChemWindow
        self.windows.append(self.structure_viewer)  # Add the ChemWindow to our list of open windows

if __name__ == '__main__':
    app = QApplication(sys.argv)  # Create a new instance of the QApplication class
    ex = Launcher()  # Create a new instance of our Launcher class
    ex.show()  # Show the Launcher window
    sys.exit(app.exec_())  # Start the application's event loop
