# main.py
import sys
from PyQt5.QtWidgets import QApplication
from launcher import Launcher
from db_core import init_connection_pool  # Import the function to initialize the connection pool

if __name__ == '__main__':
    # Initialize the connection pool before creating the Launcher
    # Replace 1 and 10 with your desired minimum and maximum connections
    init_connection_pool(1, 10)

    app = QApplication(sys.argv)
    launcher = Launcher()
    launcher.show()
    sys.exit(app.exec_())
