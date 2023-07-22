# Alembiq - A Computational Chemistry Workbench

Welcome to the Alembiq repository! This repository serves as a work in progress collection of tools and programs for computational chemistry. The workbench is powered by RDKit, Flask, and NGL Viewer.
## Version

The current release version of the workbench is v0.5

## About

The workbench currently contains four programs:
- A molecule viewer that can render molecules from the database, or draw new ones from SMILES or InChI strings.
- A molecule inserter that allows a user to add a molecule to the database using SMILES or InChI, and computes the molblock.
- A database initializer that sets up a database with model data
  - Model data: ethanol, water, oxygen, carbon dioxide and a reaction entry describing the combustion of ethanol
- A database browser/editor

Please note that this workbench is still under development and may not have all the features and capabilities you might expect from a complete computational chemistry toolset. However, we are actively working on expanding its functionality and adding more tools in the future.

## Getting Started

### Database Schema

See: [Database Schema](./db_schema.md)


To use the workbench, you will need to have Python installed on your server. Additionally, make sure you have the required dependencies listed in the `requirements.txt` file.

1. Clone this repository to your local machine.
2. Install the required dependencies:

   ```
   pip install -r requirements.txt
   ```

3. Run the launcher:

   ```
   python main.py
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
