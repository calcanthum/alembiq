from flask import Flask, jsonify, make_response, render_template, request, redirect, url_for
from db_init import DBInitializer
from db_core import init_connection_pool, get_connection_from_pool, release_connection_to_pool, fetch_molecules, connection_pool
from rdkit import Chem
from mol_insert import insert_molecule
import uuid
from db_browser import fetch_table_names, fetch_table_data, update_db

# Import the molecule module
from moldraw import get_molecule_structure

app = Flask(__name__)

@app.route('/')
def home():
    return render_template("index.html")

@app.route('/moldraw')
def molview():
    return render_template("moldraw.html")

@app.route('/db_init')
def db_init():
    try:
        print("Starting database initialization...")
        db_initializer = DBInitializer()  # Initialize the database using DBInitializer class
        print("Database initialization completed successfully.")
        response = make_response("Database initialized successfully.")
        response.mimetype = "text/plain"  # Set the response mimetype to plain text
        return response
    except Exception as e:
        return str(e)

# DB Browser related functions and routes


@app.route('/db_browser', methods=['GET', 'POST'])
def db_browser():
    # Get a connection from the pool
    conn = get_connection_from_pool()

    # Get the selected table name from the form (if any)
    selected_table = request.form.get('selected_table')

    # Fetch available table names from the database
    table_names = fetch_table_names(conn)

    # Fetch data from the selected table or the default table 'molecules'
    if selected_table and selected_table in table_names:
        data, column_names = fetch_table_data(conn, selected_table)
    else:
        data, column_names = fetch_molecules(conn)

    # Release the connection back to the pool
    release_connection_to_pool(conn)

    # Pass the data and table names to the HTML template
    return render_template('db_browser.html', data=data, table_names=table_names, selected_table=selected_table)

@app.route('/update', methods=['POST'])
def update():
    update_db(request)
    return redirect(url_for('db_browser'))


@app.route('/molecule', methods=['GET', 'POST'])
@app.route('/molecule/<id>', methods=['GET'])
def molecule(id=None):
    if request.method == 'POST':
        molecule_string = request.json.get('molecule_string', '')
        return get_molecule_structure(molecule_string)
    elif id:
        try:
            conn = get_connection_from_pool()
            cursor = conn.cursor()
            cursor.execute("SELECT smiles FROM molecules WHERE id = %s", (id,))
            molecule_string = cursor.fetchone()[0]
            release_connection_to_pool(conn)
        except Exception as e:
            print(f"Error retrieving molecule: {e}")
            return make_response(jsonify({"error": "Error retrieving molecule."}), 500)

        return get_molecule_structure(molecule_string)
    else:
        try:
            conn = get_connection_from_pool()
            molecules = fetch_molecules(conn)
            release_connection_to_pool(conn)
            return jsonify(molecules)
        except Exception as e:
            print(f"Error retrieving molecules: {e}")
            return make_response(jsonify({"error": "Error retrieving molecules."}), 500)

@app.route('/molecule_identifiers', methods=['GET'])
def get_molecule_identifiers():
    try:
        conn = get_connection_from_pool()
        cursor = conn.cursor()
        cursor.execute("SELECT id, smiles FROM molecules")
        identifiers = cursor.fetchall()
        release_connection_to_pool(conn)
        return jsonify(identifiers)
    except Exception as e:
        print(f"Error retrieving molecule identifiers: {e}")
        return make_response(jsonify({"error": "Error retrieving molecule identifiers."}), 500)



@app.route('/molecule_insert', methods=['POST'])
def molecule_insert():
    return insert_molecule(request)

if __name__ == '__main__':
    try:
        # Initialize the connection pool with appropriate minconn and maxconn values
        init_connection_pool(minconn=1, maxconn=10)

        # Run the app with debug mode
        app.run(debug=True)
    finally:
        # Close all connections when the app is shutting down
        if connection_pool:
            connection_pool.closeall()