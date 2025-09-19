# Parse any chemical file

This repo is built in order to provide a python file to parse most of molecular files.

Format supported : 

- .mol – MDL Molfile
- .sdf – Structure Data File (molecule sets, still small molecules)
- .mol2 – Tripos Mol2
- .cml – Chemical Markup Language
- .xyz – Cartesian coordinates
- .mae – Maestro small-molecule format
- .sd – SDfile variant
- .gjf / .com – Gaussian input (geometry of one molecule)
- .chk – Gaussian checkpoint (molecular state data)
- .fchk – Gaussian formatted checkpoint
- .cube – Gaussian volumetric data tied to a single molecule
- .mop – MOPAC input
- .mopout – MOPAC output
- .pdbqt – AutoDock input ligand format (derived from PDB, for small molecules only)
- .pqr – PDB with charges and radii (often for small molecules)
- .inchi – IUPAC InChI string saved as file
- .inchikey – Hashed InChI representation
- .smiles / .smi – SMILES line notation file
- .json (chem-specific schemas like ChemJSON or PubChem JSON for molecules)
- .yaml (used in some molecular toolkits for small-molecule storage)
- .mrv – ChemAxon Marvin format
- .cdx / .cdxml – ChemDraw binary and XML molecule formats
- .rxn

3 Steps will be used to return the final object :

## Detection of the extension

Either by a rule based (RB) or a simple Machine Learning (ML) model. This step is still to be decided.

## Parsing the data

RB or ML will apply the different parsing function, which will return those numpy array : 4*number_of_atoms

If it is not a text based represensation like smiles or inchi

--------------------------------------
| Atom | X_Coord | Y_Coord | Z_Coord | 
--------------------------------------

Each will parsing function will be coupled with a lenght detector to check all molecules' atom have the expected lenght

## Conversion to rdkit object

RDKIT is one of the most used open sourced chemical package avaialble, all molecule will final be converted to a rdkit mol object.
