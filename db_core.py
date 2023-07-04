import os
from dotenv import load_dotenv
from psycopg2 import connect, OperationalError, DatabaseError, sql
from psycopg2.pool import SimpleConnectionPool

load_dotenv()

DB_NAME = os.getenv('DB_NAME')
DB_USER = os.getenv('DB_USER')
DB_PASSWORD = os.getenv('DB_PASSWORD')
DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')

connection_pool = None

def init_connection_pool(minconn, maxconn):
    global connection_pool
    if not connection_pool:
        connection_pool = SimpleConnectionPool(minconn, maxconn, database=DB_NAME, user=DB_USER, password=DB_PASSWORD, host=DB_HOST, port=DB_PORT)
        if connection_pool:
            print("Successfully created connection pool")
        else:
            print("Failed to create connection pool")

def get_connection_from_pool():
    if connection_pool is None:
        raise Exception("Connection pool is not initialized")

    return connection_pool.getconn()

def release_connection_to_pool(conn):
    return connection_pool.putconn(conn)

def fetch_molecules(conn):
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT smiles, iupac FROM molecules")
        rows = cursor.fetchall()
        return rows
    except DatabaseError as e:
        print(f"Database error: {e}")
        conn.close()
        return []

def add_molecule(conn, molecule_type, other_type, molecule_string, other_string, molblock, friendly_name=None,
                 iupac_name=None, molecule_uuid=None):
    try:
        cursor = conn.cursor()

        values = (molecule_string, other_string, molblock, friendly_name, iupac_name, molecule_uuid)

        if friendly_name is None and iupac_name is None:
            insert_sql = f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, uuid) VALUES (%s, %s, %s, %s)"
            cursor.execute(insert_sql, values[:-2])  # Exclude the last two elements from values tuple
        elif friendly_name is not None and iupac_name is None:
            insert_sql = f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, friendly, uuid) VALUES (%s, %s, %s, %s, %s)"
            cursor.execute(insert_sql, (*values[:-1],))  # Exclude the last element from values tuple
        elif friendly_name is None and iupac_name is not None:
            insert_sql = f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, iupac, uuid) VALUES (%s, %s, %s, %s, %s)"
            cursor.execute(insert_sql, (*values[:-1],))  # Exclude the last element from values tuple
        else:
            insert_sql = f"INSERT INTO molecules ({molecule_type}, {other_type}, molblock, friendly, iupac, uuid) VALUES (%s, %s, %s, %s, %s, %s)"
            cursor.execute(insert_sql, values)

        conn.commit()
        cursor.close()
        return True
    except DatabaseError as e:
        print(f"Database error: {e}")
        conn.rollback()
        return False
