from flask import Flask, render_template, request, redirect, url_for
from db_core import init_connection_pool, get_connection_from_pool, release_connection_to_pool, fetch_molecules


app = Flask(__name__)

# Set your desired min and max connections for the connection pool
MIN_CONNECTIONS = 1
MAX_CONNECTIONS = 10

# Initialize the connection pool
init_connection_pool(MIN_CONNECTIONS, MAX_CONNECTIONS)


@app.route('/', methods=['GET', 'POST'])
def index():
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

def fetch_table_names(conn):
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='public';")
        rows = cursor.fetchall()
        return [row[0] for row in rows]
    except Exception as e:
        print(f"Error fetching table names: {e}")
        return []

def fetch_table_data(conn, table_name):
    try:
        cursor = conn.cursor()
        # Added ORDER BY clause here to always order rows by id
        cursor.execute(f"SELECT * FROM {table_name} ORDER BY id")
        rows = cursor.fetchall()
        column_names = [desc[0] for desc in cursor.description]
        data = [dict(zip(column_names, row)) for row in rows] if rows else []
        return data, column_names
    except Exception as e:
        print(f"Error fetching data from table '{table_name}': {e}")
        return [], []  # Return empty lists in case of error

def update_db(request):
    conn = get_connection_from_pool()
    cursor = conn.cursor()
    cursor.execute(
        f"UPDATE {request.form.get('table')} SET {request.form.get('column')} = %s WHERE id = %s",
        (request.form.get('value'), request.form.get('id'))
    )
    conn.commit()
    release_connection_to_pool(conn)

@app.route('/update', methods=['POST'])
def update():
    conn = get_connection_from_pool()
    cursor = conn.cursor()
    cursor.execute(
        f"UPDATE {request.form.get('table')} SET {request.form.get('column')} = %s WHERE id = %s",
        (request.form.get('value'), request.form.get('id'))
    )
    conn.commit()
    release_connection_to_pool(conn)
    return redirect(url_for('index'))


if __name__ == '__main__':
    app.run(debug=True)
