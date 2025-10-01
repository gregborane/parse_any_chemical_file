# Parse any chemical file

## dependencies

```
pip install numpy pyyaml json
```

This repo is built in order to provide a python file to parse most of molecular files.

Format supported : 

| Format        | Description                                                               | Is done |
| --------------|---------------------------------------------------------------------------|----------
|.mol           | MDL Molfile                                                               | Y(bonds)|
|.sdf           | Structure Data File (molecule sets, still small molecules)                | Y       |
|.mol2          | Tripos Mol2                                                               | Y       |
|.cml           | Chemical Markup Language                                                  | Y       |
|.xyz           | Cartesian coordinates                                                     | Y       |
|.mae           | Maestro small-molecule format                                             | N       |
|.sd            | SDfile variant                                                            | N       |
|.gjf / .com    | Gaussian input (geometry of one molecule)                                 | N       |
|.chk           | Gaussian checkpoint (molecular state data)                                | N       |
|.fchk          | Gaussian formatted checkpoint                                             | N       |
|.cube          | Gaussian volumetric data tied to a single molecule                        | N       |
|.mop           | MOPAC input                                                               | N       |
|.mopout        | MOPAC output                                                              | N       |
|.pdbqt         | AutoDock input ligand format (derived from PDB, for small molecules only) | N       |
|.pqr           | PDB with charges and radii (often for small molecules)                    | N       |
|.inchi         | IUPAC InChI string saved as file                                          | N       |
|.inchikey      | Hashed InChI representation                                               | N       |
|.smiles / .smi | SMILES line notation file                                                 | N       |
|.json          | chem-specific schemas like ChemJSON or PubChem JSON for molecules)        | Y       |
|.yaml          | (used in some molecular toolkits for small-molecule storage)              | Y       |
|.mrv           | ChemAxon Marvin format                                                    | N       |
|.cdx / .cdxml  | ChemDraw binary and XML molecule formats                                  | N       |
|.rxn           | can contain one molecule in some case                                     | N       |

3 Steps will be used to return the final object :

## Detection of the extension

First using a rule based (RB).\
NOT YET Later using a simple Machine Learning (ML) model to speed up proccess and code readability.

## Parsing the data

RB or ML will apply the different parsing function, which will return those numpy array : 4*number_of_atoms

If it is not a text based represensation like smiles or inchi

| Atom | X_Coord | Y_Coord | Z_Coord | 
|------|---------|---------|---------|

Each will parsing function will be coupled with a lenght detector to check all molecules' atom have the expected lenght

## Conversion to rdkit object

RDKIT is one of the most used open sourced chemical package avaialble, all molecule will final be converted to a rdkit mol object for easier post treatment.

## Final output

A dict with : 

{\
"Bonds"       : (if present take indices in file, then determined),\
"Double"      : (if present take indices in file, then determined),\
"Triple"      : (if present take indices in file, then determined),\
"Properties"  : dict(properties found in the file if any),\
"Coordinates" : np.ndarray(atoms, coordX, coordY, coordZ),\
"RDKITMol"    : RDKIT.mol\
}
