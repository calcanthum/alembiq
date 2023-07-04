import os
import psycopg2
from dotenv import load_dotenv
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QTableWidget, QTableWidgetItem, QPushButton
from db_core import get_connection_from_pool, release_connection_to_pool

load_dotenv()

DB_NAME = os.getenv("DB_NAME")
DB_USER = os.getenv("DB_USER")
DB_PASSWORD = os.getenv("DB_PASSWORD")
DB_HOST = os.getenv("DB_HOST")
DB_PORT = os.getenv("DB_PORT")

class DBBrowser(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self.table = QTableWidget()
        self.delete_button = QPushButton("Delete Selected")
        self.save_button = QPushButton("Save Changes")

        # Save button is disabled by default
        self.save_button.setEnabled(False)

        # Rows marked for deletion
        self.rows_to_delete = {}

        self.layout.addWidget(self.table)
        self.layout.addWidget(self.delete_button)
        self.layout.addWidget(self.save_button)
        self.setLayout(self.layout)

        self.load_data()

        self.delete_button.clicked.connect(self.delete_row)
        self.save_button.clicked.connect(self.save_changes)

    def load_data(self):
        conn = get_connection_from_pool()
        cur = conn.cursor()
        cur.execute('SELECT id, smiles, molblock, inchi, friendly, iupac FROM molecules')
        data = cur.fetchall()
        cur.close()
        release_connection_to_pool(conn)

        self.table.setRowCount(0)
        self.table.setColumnCount(6)  # Update column count to include 'friendly' and 'iupac'
        self.table.setHorizontalHeaderLabels(['ID', 'Smiles', 'MolBlock', 'InChI', 'Friendly', 'IUPAC'])  # Update header labels
        for row_data in data:
            row_number = self.table.rowCount()
            self.table.insertRow(row_number)
            for column_number, data in enumerate(row_data):
                self.table.setItem(row_number, column_number, QTableWidgetItem(str(data)))

    def delete_row(self):
        current_row = self.table.currentRow()
        cell_item = self.table.item(current_row, 0)
        if cell_item is not None:
            # Add row ID to rows_to_delete dictionary
            self.rows_to_delete[current_row] = int(cell_item.text())
            # Remove row from view
            self.table.removeRow(current_row)
            # Enable save button
            self.save_button.setEnabled(True)

    def save_changes(self):
        conn = get_connection_from_pool()
        cur = conn.cursor()
        for row_id in self.rows_to_delete.values():
            cur.execute('DELETE FROM molecules WHERE id = %s', (row_id,))
        conn.commit()
        cur.close()
        release_connection_to_pool(conn)
        # Clear rows_to_delete and reload data
        self.rows_to_delete.clear()
        self.load_data()
        # Disable save button
        self.save_button.setEnabled(False)

if __name__ == "__main__":
    app = QApplication([])
    db_browser = DBBrowser()
    db_browser.show()
    app.exec_()
