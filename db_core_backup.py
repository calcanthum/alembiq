import os
from dotenv import load_dotenv
from psycopg2 import DatabaseError
from psycopg2.pool import SimpleConnectionPool
import traceback

load_dotenv()

DB_NAME = os.getenv('DB_NAME')
DB_USER = os.getenv('DB_USER')
DB_PASSWORD = os.getenv('DB_PASSWORD')
DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')

connection_pool = None

def execute_sql(conn, sql, params):
    try:
        with conn.cursor() as cursor:
            cursor.execute(sql, params)
            conn.commit()
    except Exception as e:
        print(f"Error executing SQL statement: {e}")
        conn.rollback()

def execute_sql_with_error_handling(conn, sql, params):
    try:
        execute_sql(conn, sql, params)
    except Exception as e:
        error_info = traceback.format_exc()  # Get the full traceback
        print(f"Error executing SQL statement: {e}\n{error_info}")
        raise e

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
