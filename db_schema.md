# Alembiq Database Schema

The Alembiq Database is designed to store and organize chemical data related to molecules, reactions, and catalysts. The schema consists of the following tables:

---
## `molecules`
- `id`: An identifier for each molecule.
- `uuid`: A universally unique identifier for each molecule.
- `smiles`: The Simplified Molecular Input Line Entry System (SMILES) notation for the molecule.
- `inchi`: The International Chemical Identifier (InChI) for the molecule.
- `iupac`: The International Union of Pure and Applied Chemistry (IUPAC) name of the molecule.
- `smarts`: SMILES arbitrary target specification (SMARTS) string of the molecule
---
## `molecule_properties`
- `id`: An identifier for each molecule property entry.
- `molecule_id`: A foreign key linking to the corresponding molecule's `id` in the `molecules` table.
- `molecular_weight`: The molecular weight of the molecule.
- `rotatable_bond_count`: The count of rotatable bonds in the molecule.
- `aromaticity`: Information about whether the molecule is aromatic or not.
---
## `molblocks`
- `id`: An identifier for each molblock entry.
- `molecule_id`: A foreign key linking to the corresponding molecule's `id` in the `molecules` table.
- `molblock`: The molecular structure represented in the 'molblock' format.
---
## `reaction_participants`
- `id`: An identifier for each reaction participant entry.
- `reaction_id`: A foreign key linking to the corresponding reaction's `id` in the `reactions` table.
- `molecule_id`: A foreign key linking to the corresponding molecule's `id` in the `molecules` table.
- `catalyst_id`: A nullable foreign key linking to the corresponding catalyst's `id` in the `catalysts` table, if the molecule acts as a catalyst.
- `stoichiometry`: The stoichiometry of the molecule in the reaction.
- `role`: The role of the molecule in the reaction (e.g., 'reactant', 'product', 'solvent', 'catalyst').
---
## `reaction_conditions`
- `id`: An identifier for each reaction condition entry.
- `reaction_id`: A foreign key linking to the corresponding reaction's `id` in the `reactions` table.
- `temp_kelvin`: The temperature of the reaction in Kelvin.
- `pressure_kpa`: The pressure of the reaction in kilopascals.
- `ph`: The pH level of the reaction.
- `duration_seconds`: The duration of the reaction in seconds.
- `reaction_mechanism`: Information about the reaction mechanism.
---
## `catalysts`
- `id`: An identifier for each catalyst entry.
- `molecule_id`: A nullable foreign key linking to the corresponding molecule's `id` in the `molecules` table, if the catalyst is a molecule.
- `name`: The name of the catalyst.
- `description`: A text field to describe the catalyst (nullable).
- `catalyst_type`: The type of catalyst.
- `catalyst_details`: Additional details about the catalyst.
---
## `reactions`
- `id`: An identifier for each reaction.
- `uuid`: A universally unique identifier for each reaction.
- `reaction_type`: The classification or type of the reaction.
- `rxn`: RXN file for the reaction.