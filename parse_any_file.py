# %%
import numpy as np
import json
import yaml
import xml.etree.ElementTree as ET
from rdkit import Chem

def check_len(molecule: dict) -> bool:
    
    if molecule is not None:
        return True

def test_extensions(content: str|list, extension: str) -> bool:
        if isinstance(content, str) and content.split(".")[-1].lower() != extension.lower():
            raise ValueError("File is not " + extension)
        else :
            return True

def read_lines(content : str, extension: str) -> list:
    test_extensions(content, extension)
    if isinstance(content, str):
        with open(content, "r") as f:
            return f.readlines()

# %%
class XYZ:
    """
    accept as valid input :

        - direct path to the file 
        - ""readlines"" or list if same building is used
    
    """

    @staticmethod
    def get_num_atoms(content : str|list) -> int:
        
        lines = read_lines(content, "xyz") if isinstance(content, str) else content
        try:
            return int(lines[0])
        except (IndexError, ValueError):
            raise ValueError("Fail parsing num atom: error line 1")

    @staticmethod
    def get_comment(content : str|list) -> str:
        
        lines = read_lines(content, "xyz") if isinstance(content, str) else content
        try:
            return lines[1].strip()
        except IndexError:
            raise ValueError("Fail parsing comment: error line 2")

    @staticmethod
    def get_coord(content: str|list) -> tuple:

        lines = read_lines(content, "xyz") if isinstance(content, str) else content
        try:
            a = [line.split()[0] for line in lines[2:]]
            x = [line.split()[1] for line in lines[2:]]
            y = [line.split()[2] for line in lines[2:]]
            z = [line.split()[3] for line in lines[2:]]
            return a, np.array([x, y, z], dtype=np.float64)
        except (IndexError, ValueError):
            raise ValueError(f"Fail parsing coord: error within line 3 to {len(lines)+1}")

    @staticmethod
    def parse_xyz(content: str|list) -> dict:

        lines = read_lines(content, "xyz") if isinstance(content, str) else content
        return {
            "num_atom" : XYZ.get_num_atoms(lines),
            "comment"  : XYZ.get_comment(lines),
            "atoms"    : XYZ.get_coord(lines)[0],
            "coords"   : XYZ.get_coord(lines)[1],
        }

# %%
class MOL:
    """
    Accept as valid input:
        - Path to file
        - List from "readlines command"

    """
    
    @staticmethod
    def get_num_info(content: str|list) -> tuple:
         
        lines = read_lines(content, "mol") if isinstance(content, str) else content
        
        try:
            if "V3000" in lines[3]:
                raise ValueError("V3000 formatting is not supported")
                
            counts = lines[3].split()
            num_atom = int(counts[0])
            num_bonds = int(counts[1])
            return num_atom, num_bonds
        except (ValueError, IndexError):
            raise ValueError("Fail parsing num bonds or atoms at line 4")

    @staticmethod
    def get_coord(content: str|list) -> tuple:
        
        lines = read_lines(content, "mol") if isinstance(content, str) else content
        num_atom = MOL.get_num_info(lines)[0]
        
        try:
            a = [line.split()[3] for line in lines[4:4+num_atom]]
            x = [line.split()[0] for line in lines[4:4+num_atom]]
            y = [line.split()[1] for line in lines[4:4+num_atom]]
            z = [line.split()[2] for line in lines[4:4+num_atom]]
            return a, np.array([x,y,z], dtype=np.float64) 
        except (ValueError, IndexError):
            raise ValueError(f"Error parsing coords between lines 5;{num_atom+1}")

    @staticmethod
    def get_bonds(content: str|list) -> tuple:

        lines = read_lines(content, "mol") if isinstance(content, str) else content
        num_atom, num_bonds = MOL.get_num_info(content) 
        
        try:
            bonds, double_bonds, triple_bonds = [], [], []
            for line_bond in lines[4+num_atom:4+num_atom+num_bonds]:
                tokens = line_bond.split()
                a1, a2, bond_type = int(tokens[0])-1, int(tokens[1])-1, int(tokens[2])-1
                if bond_type == 1:
                    bonds.append(
                        (a1, a2)
                    )
                elif bond_type == 2:
                    double_bonds.append(
                        (a1, a2)
                    )
                elif bond_type == 3:
                    triple_bonds.append(
                        (a1, a2)
                    )
            return bonds, double_bonds, triple_bonds
        except (ValueError, IndexError):
            raise ValueError(f"Error parsing bonds information between lines {5+num_atom} and {5+num_atom+num_bonds}")

    @staticmethod
    def get_information(content: str|list) -> dict:
    
        lines = read_lines(content, "mol") if isinstance(content, str) else content
        num_atom, num_bonds = MOL.get_num_info(content)
        
        try:
            properties = {}
            for line in lines[4+num_atom+num_bonds:]:
                if line.startswith("M  END"):
                    break
                tokens = line.split()
                if len(tokens) > 2:
                    properties[tokens[1]] = tokens[2:]
            return properties
        except (ValueError, IndexError):
            raise ValueError(f"Error Parsing Information between lines {6+num_atom+num_bonds} and {len(lines)+1}")

    @staticmethod
    def parse_mol(content: str|list) -> dict:
    
        lines = read_lines(content, "mol") if isinstance(content, str) else content
        
        return {
            "num_atom"     : MOL.get_num_info(lines)[0],
            "properties"   : MOL.get_information(lines),
            "bonds"        : MOL.get_bonds(lines)[0],
            "double_bonds" : MOL.get_bonds(lines)[1],
            "triple_bonds" : MOL.get_bonds(lines)[2],
            "atoms"        : MOL.get_coord(lines)[0],
            "coord"        : MOL.get_coord(lines)[1],
        }

# %%
def parse_mol2(path: str) -> dict:
    """
    Parse MOL2 like XYZ adding bonds
    information since they exist.
    """

    with open(path, "r") as mol2:
        lines = mol2.readlines()

    bonds, double_bonds, triple_bonds = [], [], []
    coords, charges = [], []
    mol_type, comment, charges_type = None, None, None
    num_atom, num_bond = 0, 0

    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>MOLECULE"):
            num_atom, num_bond = map(int, lines[i + 2].split()[:2])
            mol_type = lines[i + 3].strip()
            charges_type = lines[i + 4].strip()
            if len(lines) > i + 5 and not lines[i+5].startswith("@<TRIPOS>"):
                comment = lines[i + 5].strip()        

        elif line.startswith("@<TRIPOS>ATOM"):
            for line_coord in lines[i + 1 : i + 1 + num_atom]:
                tokens = line_coord.split()
                atoms = tokens[1]
                x, y, z = map(float, tokens[2:5])
                coords.append((atoms, x, y, z))
                charges.append(float(tokens[-1]))

        elif line.startswith("@<TRIPOS>BOND"):
            for line_bond in lines[i + 1 : i + 1 + num_bond]:
                tokens = line_bond.split()
                a1, a2, bond_type = tokens[1], tokens[2], tokens[3]
                if bond_type == "1":
                    bonds.append((a1, a2))
                elif bond_type in ("2", "ar"):
                    double_bonds.append((a1, a2))
                elif bond_type == "3":
                    triple_bonds.append((a1, a2))
        else:
            raise ValueError(
                "Cannot parse chemical information without \n MOLECULE, ATOM, BOND sections."
            )

    molecule = {
        "num_atom": num_atom,
        "num_bond": num_bond,
        "mol_type": mol_type,
        "charges_type": charges_type,
        "comment": comment,
        "coords": np.array(coords),
        "charges": charges,
        "bonds": bonds,
        "double_bonds": double_bonds,
        "triple_bonds": triple_bonds,
    }
    return molecule


def parse_json(path: str):
    """
    JSON file do not have systematic representation.
    So they are just parsed using built in python parser.
    """

    with open(path, "r") as file:
        data = json.load(file)

    return data


def parse_yaml(path: str):
    """
    YAML file do not have systematic representation.
    So they are just parsed using built in function.
    """

    with open(path, "r") as file:
        data = yaml.safe_load_all(file)

    return data

def parse_sdf(path: str) -> list:
    """
    Parse SDF file into list of molecule dictionaries.
    Each molecule contains atoms, bonds, and property data.
    """

    with open(path, "r") as sdf:
        content = sdf.read().strip()

    molecules = []
    blocks = content.split("$$$$")
    for block in blocks:
        lines = [l for l in block.strip().splitlines() if l.strip()]
        if not lines:
            continue

        # detect V3000 unsupported
        if any("V3000" in l for l in lines):
            raise KeyError("V3000 formatting is not supported")

        # header: 3 lines + counts line
        counts_line = lines[3]
        num_atom = int(counts_line.split("")[0])
        num_bond = int(counts_line.split()[1])

        atom_start = 4
        atom_end = atom_start + num_atom
        bond_end = atom_end + num_bond

        atoms = []
        for l in lines[atom_start:atom_end]:
            parts = l.split()
            x, y, z = map(float, parts[0:3])
            atom_symbol = parts[3]
            atoms.append((atom_symbol, x, y, z))

        bonds, double_bonds, triple_bonds = [], [], []
        for l in lines[atom_end:bond_end]:
            parts = l.split()
            a1, a2, btype = parts[0], parts[1], parts[2]
            if btype == "1":
                bonds.append((a1, a2))
            elif btype in ("2", "ar"):
                double_bonds.append((a1, a2))
            elif btype == "3":
                triple_bonds.append((a1, a2))

        properties = {}
        i = bond_end
        while i < len(lines):
            if lines[i].startswith(">"):
                key = lines[i].strip()[3:-1]
                val = []
                i += 1
                while i < len(lines) and not lines[i].startswith(">"):
                    val.append(lines[i].strip())
                    i += 1
                properties[key] = "\n".join(val)
            else:
                i += 1

        molecule = {
            "num_atom": num_atom,
            "num_bond": num_bond,
            "coord": np.array(atoms),
            "bonds": bonds,
            "double_bonds": double_bonds,
            "triple_bonds": triple_bonds,
            "properties": properties,
        }
        molecules.append(molecule)

    return molecules

def parse_cml(path: str) -> dict:
    """
    Parse cml file which are a subclass of xml files.

    In the doc, there are much more formatting possible so I only picked the most common ones :
    - Atoms
    - Coordinates
    - Bonds Information
    - Properties

    To add more section parsing you can try:
    for x in root.findall(.//cml:{section}):
        x.get("attributes")

    And associated code to parse data section.

    More over, this parsing is built over cml.org schematic.

    """
    ns = {"cml": "http://www.xml-cml.org/schema"}
    tree = ET.parse(path)
    root = tree.getroot()

    atoms = [atom.get("elementType") for atom in root.findall(".//cml:atom", ns)]

    # Element's coordinate can be given in either two dimensions or three dimensions.
    try:
        x = [atom.get("x3") for atom in root.findall(".//cml:atom", ns)]
        y = [atom.get("y3") for atom in root.findall(".//cml:atom", ns)]
        z = [atom.get("z3") for atom in root.findall(".//cml:atom", ns)]

    finally:
        x = [atom.get("x2") for atom in root.findall(".//cml:atom", ns)]
        y = [atom.get("y2") for atom in root.findall(".//cml:atom", ns)]

    premol = np.array([atoms, x, y, z])
    num_atom = len(atoms)

    bonds, double_bonds, triple_bonds = list(), list(), list()

    for atom_bonds in root.findall(".//cml:bonds", ns):
        for key in ["atomRefs", "atomRefs2", "atomRefs3", "atomRefs4"]:

            # atomRefs# can be used differentyly depending on the atom si it is
            # a necessity to test them all, might lead to an error
            try:
                if atom_bonds.get("order") == "1" and atom_bonds.get(key):
                    bond = str(atom_bonds.get(key))
                    bond_atom = bond.split(" ")
                    parsed = [int(bond_at.lstrip("a")) - 1 for bond_at in bond_atom]
                    bonds.append(tuple(parsed))

                elif atom_bonds.get("order") == "2" and atom_bonds.get(key):
                    bond = str(atom_bonds.get(key))
                    bond_atom = bond.split(" ")
                    parsed = [int(bond_at.lstrip("a")) - 1 for bond_at in bond_atom]
                    double_bonds.append(tuple(parsed))

                elif atom_bonds.get("order") == "3" and atom_bonds.get(key):
                    bond = str(atom_bonds.get(key))
                    bond_atom = bond.split(" ")
                    parsed = [int(bond_at.lstrip("a")) - 1 for bond_at in bond_atom]
                    triple_bonds.append(tuple(parsed))

            finally:
                pass

    properties = list()
    for prop in root.findall(".//cml:property", ns):
        title = prop.get("title")
        scalar = prop.find("cml:scalar", ns)
        if scalar is not None:
            value = scalar.text
            property = {title: value}
        else:
            continue
        properties.append(property)

    molecule = {
        "num_atom": num_atom,
        "properties": properties,
        "coord": premol,
        "bonds": bonds,
        "double_bonds": double_bonds,
        "triple_bonds": triple_bonds,
    }

    return molecule

def parse_smi(path : str) -> dict:
    """
    Parse smiles files and convert them into RDKIT mol object
    """

    with open(path, "r") as smi:
        line = smi.readlines()[0]
    
    molobject = Chem.MolFromSmiles(line)
    molecule = {
            "molecule"  : line,
            "MolObject" : molobject
    }

    return molecule

def parse_inchi(path : str) -> dict:
    """
    Parse inchi file and convert them into RDKIT mol object
    """

    with open(path, "r") as inchi:
        line = inchi.readlines()[0]

    molobject = Chem.inchi.MolFromInchi(line)
    molecule = {
            "text"      : line,
            "molobject" : molobject
    }

    return molecule

"""
class ReadGaussian:

    @static
    def read_lines():
        
        Get Data
       
        if test_extensions(path, "log"): 
            with open(path,"r") as gaussian_file:
                    lines = gaussian_file.readlines()
                    lines = lines
        return lines

    def extract_coordinates():
      
        Obtain Atoms coordinates and atomic number
     
        lines = read_lines()
        num_atom, geometry_indice = None, None

        for i, line in enumerate(lines[::-1]):
            split_line = line.split()

            if "NAtoms=" in split_line:
                num_atom = int(line.split()[1])
            
            if "standard" in split_line and "orientation" in split_line:
                geometry_indice = i
        
        if geometry_indice is None :
            raise ValueError("No Coordinates Found")

        if num_atom is None:
            raise ValueError("Number of Atom not Found")

        start_coord = len(lines) - geometry_indice + 5
        end_coord = start_coord + num_atom

        x = [float(line.split()[3]) for line in lines[start_coord: end_coord]]
        y = [float(line.split()[4]) for line in lines[start_coord: end_coord]]
        z = [float(line.split()[5]) for line in lines[start_coord: end_coord]]
        atom = [float(line.split()[1]) for line in lines[start_coord: end_coord]]
   
        premol = np.array([atom, x, y, z])
    
        return premol

    def extract_charges(path: str):
       
        Obtain partial charges from NBO calculations
        Mulliken, 
      

        lines = read_lines(path, "log")
        indice_mulliken = None, 

        for i, lines in enumerate(lines[::-1]):
            split_line = line.split()

            if "Mulliken" in lines.split() and len(lines.split()) == 2:
                indice_mulliken = i

            if "APT" in lines.split() and len(lines.split()) == 2:

        return "bite"

    def extract_bond_order(path: str):
     
        Obtain bond order from NBO calculation
    
        return 'bite bite'

"""
