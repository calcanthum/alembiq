# AlembIQ - A Computational Chemistry Workbench

Welcome to the AlembIQ repository! This repository serves as a work in progress collection of tools and programs for computational chemistry. The workbench is powered by RDKit, Flask, and NGL Viewer.

## Version

The current release version of the workbench is v0.5.1

## About

The workbench currently contains four programs:
- A molecule viewer that can render molecules from the database, or draw new ones from SMILES or InChI strings.
- A molecule inserter that allows a user to add a molecule to the database using SMILES or InChI, and computes the molblock and SMARTS.
- A database initializer that sets up a database with model data
  - Model data: ethanol, water, oxygen, carbon dioxide and a reaction entry describing the combustion of ethanol
- A database browser/editor
  
⚠️ Warning ⚠️
---
Do not use this in production.
This code uses a flask server in debug mode by default. It also has naïve connection pooling and little to no input sanitization or parameterization. This workbench is still under development and may not have all the features and capabilities you might expect from a complete computational chemistry toolset.

## Getting Started

To use the workbench, you will need to have Python installed and the required dependencies listed in the `requirements.txt` file. A Postgresql server and database must already exist. Database connection details are read from a(n) .env file with the following structure:
```
# .env
DB_HOST='***'
DB_NAME='***'
DB_USER='***'
DB_PASSWORD='***'
DB_PORT='***'
```
The Database Initializer app will setup the database tables for you and insert model data according to this [Database Schema](./db_schema.md)
1. Clone this repository to a local directory and enter it:

    ```
    git clone https://github.com/calcanthum/alembiq
    cd alembiq
    ```

2. Install the required dependencies:

   ```
   pip install -r requirements.txt
   ```

3. Run the server:

   ```
   python alembiq_server.py
   ```

4. Open a browser window and go to:

    ```
    127.0.0.1:5000
    ```

## Contributing

WIP

## Licensing

This project is licensed under the terms of the GNU Affero General Public License v3.0 (AGPLv3).

The GNU Affero General Public License is a free, copyleft license for software and other kinds of works, specifically designed to ensure cooperation with the community in the case of network server software. The licenses for most software and other practical works are designed to take away your freedom to share and change the works. By contrast, our General Public Licenses are intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users.

When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are designed to make sure that you have the freedom to distribute copies of free software (and charge for them if you wish), that you receive source code or can get it if you want it, that you can change the software or use pieces of it in new free programs, and that you know you can do these things.

To protect your rights, we need to prevent others from denying you these rights or asking you to surrender the rights. Therefore, you have certain responsibilities if you distribute copies of the software, or if you modify it: responsibilities to respect the freedom of others.

For more details, see the [LICENSE](./LICENSE) file in this repository, or [the full text of the license](https://www.gnu.org/licenses/agpl-3.0.en.html).


## Disclaimer

This workbench is provided as-is and without any warranty. Use it at your own risk. The authors and contributors of this repository shall not be liable for any damages or consequences arising from the use of this workbench.
