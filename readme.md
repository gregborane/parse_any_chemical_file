# Parse any chemical file

## Reasons

I am tired of parsing manually files after each project so here it is a parser for file type.

Long term project is to contain all extensions with the bare minimum extractions informations : 

- Atoms,  
- Coordinates,
- Bonds,
- Charges

Containing :

- Auto detection of extensions
- Possibility to choose which information to extract from the file (default all) 

## dependencies

```
pip install numpy pyyaml json rdkit
```

This repo is built in order to provide a python file to parse most of molecular files.

Format supported : 

| Format              | Description                                                               | Is done |
|---------------------|---------------------------------------------------------------------------|----------
|.mol                 | MDL Molfile                                                               | Y       |
|.sdf                 | Structure Data File (molecule sets, still small molecules)                | Y       |
|.mol2                | Tripos Mol2                                                               | Y       |
|.cml                 | Chemical Markup Language                                                  | Y       |
|.xyz                 | Cartesian coordinates                                                     | Y       |
|.mae                 | Maestro small-molecule format                                             | N       |
|.sd                  | SDfile variant                                                            | N       |
|.gjf / .com/ .out    | Gaussian input (geometry of one molecule)                                 | ~       |
|.cube                | Gaussian volumetric data tied to a single molecule                        | N       |
|.mop                 | MOPAC input                                                               | N       |
|.mopout              | MOPAC output                                                              | N       |
|.pdbqt               | AutoDock input ligand format (derived from PDB, for small molecules only) | N       |
|.pqr                 | PDB with charges and radii (often for small molecules)                    | N       |
|.inchi               | IUPAC InChI string saved as file                                          | Y       |
|.smiles / .smi       | SMILES line notation file                                                 | Y       |
|.json                | chem-specific schemas like ChemJSON or PubChem JSON for molecules)        | Y       |
|.yaml                | (used in some molecular toolkits for small-molecule storage)              | Y       |
|.mrv                 | ChemAxon Marvin format                                                    | N       |
|.cdx / .cdxml        | ChemDraw binary and XML molecule formats                                  | N       |

## 3 Steps will be used to return the final object :

### Detection of the extension

First using a rule based (RB).\

### Parsing the data

If it is not a text based represensation like smiles or inchi

Most of the information will be extracted and will be relative to a file.
Check docs for each extension extraction.

### Conversion to rdkit object

RDKIT is one of the most used open sourced chemical package avaialble, all molecule will final be converted to a rdkit mol object for easier post treatment.

### Final output

A dict with : 

{\
"Information Extracted",\
"RDKIT_object",\
}
