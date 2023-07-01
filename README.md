
# Alembiq - A Computational Chemistry Workbench

Welcome to the Alembiq repository! This repository serves as a work in progress collection of tools and programs for computational chemistry based on the RDKit framework and psycopg2, connecting to a PostgreSQL database (server not included).

## Version

The current release version of the workbench is v0.3.0.

## About

The workbench currently contains five programs:
- A molecule viewer that connects to a PostgreSQL database
- A database updater for adding molecules to a PostgreSQL database via SMILES strings
- A SMILES to molblock converter
- A database discrepancy checker and resolver
- A database initializer that connects to a database and creates an RDKit compatible table (`molecules`) with an ethanol molecule added.
   - The initializer will wipe any existing table called `molecules`

Please note that this workbench is still under development and may not have all the features and capabilities you might expect from a complete computational chemistry toolset. However, we are actively working on expanding its functionality and adding more tools in the future.

## Getting Started

To use the workbench, you will need to have Python installed on your system. Additionally, make sure you have the required dependencies listed in the `requirements.txt` file.

1. Clone this repository to your local machine.
2. Install the required dependencies by running the following command:

   ```
   pip install -r requirements.txt
   ```

3. Run the program using the following command:

   ```
   python main.py
   ```

   This will run the launcher application.

## Contributing

WIP

## License

WIP

## Disclaimer

This workbench is provided as-is and without any warranty. Use it at your own risk. The authors and contributors of this repository shall not be liable for any damages or consequences arising from the use of this workbench.
