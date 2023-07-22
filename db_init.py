import traceback
import uuid

import psycopg2
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

import db_core

# load environment variables
load_dotenv()
print("Starting DB Initializer...")


class DBInitializer:
    def __init__(self):
        # Get a database connection from the connection pool
        self.conn = db_core.get_connection_from_pool()
        # Create a cursor object
        self.cursor = self.conn.cursor()

        # Delete the tables if they exist, then recreate them
        self.delete_tables()
        self.create_tables()

        # Insert ethanol
        ethanol = Chem.MolFromSmiles("CCO")
        self.insert_molecule(ethanol, "Ethanol")

        # Insert water
        water = Chem.MolFromSmiles("O")
        self.insert_molecule(water, "Water")

        # Insert oxygen
        oxygen = Chem.MolFromSmiles("O=O")
        self.insert_molecule(oxygen, "Oxygen")

        # Insert carbon dioxide
        carbon_dioxide = Chem.MolFromSmiles("C(=O)=O")
        self.insert_molecule(carbon_dioxide, "Carbon dioxide")

        # Insert model reaction of the combustion of ethanol into the reactions table
        self.insert_reaction(
            reaction_type="Combustion",
            reactants=[("Ethanol", 1), ("Oxygen", 3)],
            products=[("Carbon dioxide", 2), ("Water", 3)],
        )

        # Insert reaction conditions of model reaction
        self.insert_reaction_conditions(1, 638.15, 101.325, 7, 3600, "Combustion")

        # Close the cursor and release the database connection back to the pool
        self.cursor.close()
        db_core.release_connection_to_pool(self.conn)

    def delete_tables(self):
        tables_to_drop = ["reaction_participants", "reactions", "molecule_properties", "molblocks", "molecules",
                          "reaction_conditions", "catalysts"]

        for table in tables_to_drop:
            drop_table_sql = f"DROP TABLE IF EXISTS {table} CASCADE"
            db_core.execute_sql_with_error_handling(self.conn, drop_table_sql, ())

    def create_tables(self):
        try:
            # SQL statements to create the new tables
            create_molecules_table = """
            CREATE TABLE IF NOT EXISTS molecules (
                id SERIAL PRIMARY KEY,
                uuid UUID,
                smiles TEXT,
                inchi TEXT,
                iupac TEXT
            )
            """
            create_molecule_properties_table = """
            CREATE TABLE IF NOT EXISTS molecule_properties (
                id SERIAL PRIMARY KEY,
                molecule_id INTEGER REFERENCES molecules(id),
                molecular_weight FLOAT,
                rotatable_bond_count INTEGER,
                aromaticity INTEGER
            )
            """
            create_molblocks_table = """
            CREATE TABLE IF NOT EXISTS molblocks (
                id SERIAL PRIMARY KEY,
                molecule_id INTEGER REFERENCES molecules(id),
                molblock TEXT
            )
            """
            create_reactions_table = """
            CREATE TABLE IF NOT EXISTS reactions (
                id SERIAL PRIMARY KEY,
                uuid UUID,
                reaction_type TEXT
            )
            """
            create_catalysts_table = """
            CREATE TABLE IF NOT EXISTS catalysts (
                id SERIAL PRIMARY KEY,
                molecule_id INTEGER REFERENCES molecules(id),
                name TEXT,
                description TEXT,
                catalyst_type TEXT,
                catalyst_details TEXT
            )
            """
            create_reaction_participants_table = """
            CREATE TABLE IF NOT EXISTS reaction_participants (
                id SERIAL PRIMARY KEY,
                reaction_id INTEGER REFERENCES reactions(id),
                molecule_id INTEGER REFERENCES molecules(id),
                catalyst_id INTEGER REFERENCES catalysts(id),
                stoichiometry FLOAT,
                role TEXT
            )
            """
            create_conditions_table = """
            CREATE TABLE IF NOT EXISTS reaction_conditions (
                id SERIAL PRIMARY KEY,
                reaction_id INTEGER REFERENCES reactions(id),
                temp_kelvin FLOAT,
                pressure_kpa FLOAT,
                pH FLOAT,
                duration_seconds INTEGER,
                reaction_mechanism TEXT
            )
            """

            self.cursor.execute(create_molecules_table)
            self.cursor.execute(create_molecule_properties_table)
            self.cursor.execute(create_molblocks_table)
            self.cursor.execute(create_reactions_table)
            self.cursor.execute(create_catalysts_table)  # Move this before creating "reaction_participants" table
            self.cursor.execute(create_reaction_participants_table)
            self.cursor.execute(create_conditions_table)
            self.conn.commit()
        except psycopg2.Error as e:
            error_info = traceback.format_exc()  # Get the full traceback
            print("Error Creating Tables:", e, "\n", error_info)

    def calculate_uuid(self):
        # Calculate UUID using uuid.uuid4()
        return str(uuid.uuid4())

    def insert_molecule(self, mol, iupac):
        try:
            # Generate the molecule's SMILES, molblock, InChI, and UUID
            smiles = Chem.MolToSmiles(mol)
            molblock = Chem.MolToMolBlock(mol, includeStereo=True)
            inchi = Chem.MolToInchi(mol)
            uuid = self.calculate_uuid()

            values = (uuid, smiles, inchi, iupac)
            statement = "INSERT INTO molecules (uuid, smiles, inchi, iupac) VALUES (%s, %s, %s, %s)"
            db_core.execute_sql(self.conn, statement, values)

            molecule_id = self.get_molecule_id(iupac)

            # Calculate molecule properties
            molecular_weight = Descriptors.MolWt(mol)
            rotatable_bond_count = rdMolDescriptors.CalcNumRotatableBonds(mol)
            aromaticity = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)  # returns an integer

            # Print SMILES and rotatable bond count
            print(f'SMILES: {smiles}, Rotatable Bonds: {rotatable_bond_count}')

            values = (molecule_id, molecular_weight, rotatable_bond_count, aromaticity)
            statement = "INSERT INTO molecule_properties (molecule_id, molecular_weight, rotatable_bond_count, aromaticity) VALUES (%s, %s, %s, %s)"
            db_core.execute_sql(self.conn, statement, values)

            # Insert molblock
            values = (molecule_id, molblock)
            statement = "INSERT INTO molblocks (molecule_id, molblock) VALUES (%s, %s)"
            db_core.execute_sql(self.conn, statement, values)
        except psycopg2.Error as e:
            error_info = traceback.format_exc()  # Get the full traceback
            print("Error Inserting Molecule:", e, "\n", error_info)


    def insert_reaction(self, reaction_type, reactants, products):
        try:
            # Calculate a new UUID for the reaction
            reaction_uuid = self.calculate_uuid()

            # Insert the reaction into the reactions table
            statement = "INSERT INTO reactions (uuid, reaction_type) VALUES (%s, %s)"
            self.cursor.execute(statement, (reaction_uuid, reaction_type))

            # Get the ID of the reaction we just inserted
            self.cursor.execute("SELECT id FROM reactions WHERE uuid = %s", (reaction_uuid,))
            reaction_id = self.cursor.fetchone()[0]

            # Insert the reactants and products into the reactants table
            for role, molecules in [("reactant", reactants), ("product", products)]:
                for molecule_iupac, stoichiometry in molecules:
                    # Get the molecule's ID
                    self.cursor.execute("SELECT id FROM molecules WHERE iupac = %s", (molecule_iupac,))
                    molecule_id = self.cursor.fetchone()[0]

                    # Insert the molecule as a participant in the reaction
                    statement = "INSERT INTO reaction_participants (reaction_id, molecule_id, stoichiometry, role) VALUES (%s, %s, %s, %s)"
                    self.cursor.execute(statement, (reaction_id, molecule_id, stoichiometry, role))
            self.conn.commit()
        except psycopg2.Error as e:
            error_info = traceback.format_exc()  # Get the full traceback
            print("Error Inserting Reaction:", e, "\n", error_info)

    def insert_reaction_conditions(self, reaction_id, temp_kelvin, pressure_kpa, pH, duration_seconds, reaction_mechanism):
        try:
            values = (reaction_id, temp_kelvin, pressure_kpa, pH, duration_seconds, reaction_mechanism)
            statement = "INSERT INTO reaction_conditions (reaction_id, temp_kelvin, pressure_kpa, pH, duration_seconds, reaction_mechanism) VALUES (%s, %s, %s, %s, %s, %s)"
            db_core.execute_sql(self.conn, statement, values)
        except psycopg2.Error as e:
            error_info = traceback.format_exc()  # Get the full traceback
            print("Error Inserting Reaction Conditions:", e, "\n", error_info)

    def get_molecule_id(self, iupac):
        self.cursor.execute("SELECT id FROM molecules WHERE iupac = %s", (iupac,))
        return self.cursor.fetchone()[0]


if __name__ == "__main__":
    DBInitializer()
    print("Database Initialized!")