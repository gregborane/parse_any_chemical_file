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
            raise ValueError(
            "File is not " + extension
        )
        else :
            return True

def read_lines(content : str|list, extension: str) -> list:
    test_extensions(content, extension)
    if isinstance(content, str):
        with open(content, "r") as f:
            return f.readlines()
    elif isinstance(content, list):
        return content
    else:
        raise TypeError("Input invalid type")

# %%
class XYZ:
    """
    accept as valid input :

        - direct path to the file 
        - ""readlines"" or list if same building is used
    
    """

    # Magic numbers =

    ## inside section
    a_i = 0
    x_i = 1
    y_i = 2
    z_i = 3
    
    ## sections
    comment_i  = 1
    num_atom_i = 0
    coord_i    = 2
    file_extension = "xyz"

    @staticmethod
    def get_num_atoms(content : str|list) -> int:
        
        lines = read_lines(content, XYZ.file_extension)

        try:
            return int(lines[XYZ.num_atom_i])
        except (IndexError, ValueError):
          raise ValueError(
                f"Fail parsing num atom: error line {XYZ.num_atom_i+1}"
            )

    @staticmethod
    def get_comment(content : str|list) -> str:
        
        lines = read_lines(content, "xyz") 
        try:
            return lines[XYZ.comment_i]
        except IndexError:
            raise ValueError(
                f"Fail parsing comment: error line {XYZ.comment_i+1}"
            )

    @staticmethod
    def get_coord(content: str|list) -> tuple:

        lines = read_lines(content, XYZ.file_extension)

        try:
            a = [line.split()[XYZ.a_i] for line in lines[XYZ.coord_i:]]
            x = [line.split()[XYZ.x_i] for line in lines[XYZ.coord_i:]]
            y = [line.split()[XYZ.y_i] for line in lines[XYZ.coord_i:]]
            z = [line.split()[XYZ.z_i] for line in lines[XYZ.coord_i:]]
            return a, np.array([x, y, z], dtype=np.float64)
        except (IndexError, ValueError):
            raise ValueError(
                f"Fail parsing coord: error within line {XYZ.coord_i+1} to {len(lines)+1}"
            )

    @staticmethod
    def parse_xyz(content: str|list) -> dict:

        lines = read_lines(content, XYZ.file_extension)

        return {
            "num_atom" : XYZ.get_num_atoms(lines),

            "comment" : XYZ.get_comment(lines),

            "atoms" : XYZ.get_coord(lines)[0],
            "coords" : XYZ.get_coord(lines)[1],
        }

# %%
class MOL:
    """
    Accept as valid input:
        - Path to file
        - List from "readlines command"

    """

    # Magic numbers =
    ## inside section
    a_i = 3
    x_i = 0
    y_i = 1
    z_i = 2
    atom_i = 0
    bond_i = 1

    ## sections
    num_i = 3
    file_extension = "mol"
    coord_i = 4
    
    @staticmethod
    def get_num_info(content: str|list) -> tuple:
 
        lines = read_lines(content, MOL.file_extension)

        try:
            if "V3000" in lines[MOL.num_i]:
                raise ValueError(
                "V3000 formatting is not supported"
                )
                
            counts = lines[MOL.num_i].split()
            num_atom = int(counts[MOL.atom_i])
            num_bonds = int(counts[MOL.bond_i])
            return num_atom, num_bonds
        except (ValueError, IndexError):
            raise ValueError(
                f"Fail parsing num bonds or atoms at line {MOL.num_i+1}"
            )

    @staticmethod
    def get_coord(content: str|list) -> tuple:
        
        lines = read_lines(content, MOL.file_extension)
        num_atom = MOL.get_num_info(lines)[0]
        
        try:
            a = [line.split()[MOL.a_i] for line in lines[MOL.coord_i:MOL.coord_i+num_atom]]
            x = [line.split()[MOL.x_i] for line in lines[MOL.coord_i:MOL.coord_i+num_atom]]
            y = [line.split()[MOL.y_i] for line in lines[MOL.coord_i:MOL.coord_i+num_atom]]
            z = [line.split()[MOL.z_i] for line in lines[MOL.coord_i:MOL.coord_i+num_atom]]
            return a, np.array([x,y,z], dtype=np.float64) 
        except (ValueError, IndexError):
            raise ValueError(
                f"Error parsing coords between lines {MOL.coord_i+1} and {num_atom+1}"
            )

    @staticmethod
    def get_bonds(content: str|list) -> tuple:

        lines = read_lines(content, MOL.file_extension) 
        num_atom, num_bonds = MOL.get_num_info(content) 
        
        try:
            bonds, double_bonds, triple_bonds = [], [], []
            for line_bond in lines[MOL.coord_i+num_atom:MOL.coord_i+num_atom+num_bonds]:
                tokens = line_bond.split()
                # -1 to correct index to 0
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
            raise ValueError(
                f"Error parsing bonds information between lines {MOL.coord_i+1+num_atom} and {MOL.coord_i+num_atom+num_bonds}"
            )

    @staticmethod
    def get_information(content: str|list) -> dict:
    
        lines = read_lines(content, MOL.file_extension) 
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
            raise ValueError(
                f"Error Parsing Information between lines {MOL.num_i+1+num_atom+num_bonds} and {len(lines)+1}"
            )

    @staticmethod
    def parse_mol(content: str|list) -> dict:
    
        lines = read_lines(content, MOL.file_extension) if isinstance(content, str) else content
        
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
class MOL2:
            
    # Magic numbers = 
    x_i = 2
    y_i = 3
    z_i = 4
    a_i = 1
    c_i = -1 
    
    type_i = 3
    atom_bond_num_i = 2
    comment_i = 5
    chargetypes_i = 4
    file_extension = "mol2"

    @staticmethod
    def get_num_info(content: str|list) -> tuple:

        lines = read_lines(content, MOL2.file_extension) 

        for i, line in enumerate(lines):
            if line.startswith("@<TRIPOS>MOLECULE"):
                num_atom, num_bond = map(int, lines[i + MOL2.atom_bond_num_i].split()[:MOL2.atom_bond_num_i])
                return num_atom, num_bond
        raise ValueError("No Molecule section found")

    @staticmethod
    def get_mol_type(content: str|list) -> list:

        lines = read_lines(content, MOL2.file_extension) 

        for i, line in enumerate(lines):
            if line.startswith("@<TRIPOS>MOLECULE"):
                try :
                    return lines[i + MOL2.type_i].strip()
                except (IndexError, ValueError):
                    raise ValueError(
                        f"Error Parsing molecule type at line {i+MOL2.type_i}"
                    )
        raise ValueError(
            "No MOLECULE section found"
        )

    @staticmethod
    def get_comment(content: str|list) -> list: 

        lines = read_lines(content, MOL2.file_extension) 

        for i, line in enumerate(lines):
            if len(lines) > i + MOL2.comment_i and not line.startswith("@<TRIPOS>"):
                try:
                    return lines[i + MOL2.comment_i].strip()        
                except(ValueError, IndexError):
                    raise ValueError(
                        f"Error Parsing comment at line {i+MOL2.comment_i}"
                    )
        raise ValueError(
            "No COMMENT found"
        )

    @staticmethod
    def get_charges_types(content: str|list) -> list:

        lines = read_lines(content, MOL2.file_extension) 

        for i, line in enumerate(lines):
            if line.startswith("@<TRIPOS>MOLECULE"):
                index_molecule_start = i
                try :
                    return lines[i + MOL2.chargetypes_i].strip()
                except(ValueError, IndexError):
                    raise ValueError(
                        f"Could not find charges type at line {index_molecule_start+MOL2.chargetypes_i}"
                    )
        raise ValueError(
            "No MOLECULE section found"
        )

    @staticmethod
    def get_coord(content : str|list) -> tuple:

        lines = read_lines(content, MOL2.file_extension) 
        num_atom = MOL2.get_num_info(lines)[0]

        for i, line in enumerate(lines):
            if line.startswith("@<TRIPOS>ATOM"):
                index_atom_start = i
                try:
                    # 1 skip line and 2 to reach all atoms
                    a = [line.split()[MOL2.a_i] for line in lines[i+1: num_atom]] 
                    x = [np.float64(line.split()[MOL2.x_i]) for line in lines[i+1: i+2+num_atom]]
                    y = [np.float64(line.split()[MOL2.y_i]) for line in lines[i+1: i+2+num_atom]]
                    z = [np.float64(line.split()[MOL2.z_i]) for line in lines[i+1: i+2+num_atom]]
                    c = [np.float64(line.split()[MOL2.c_i]) for line in lines[i+1: i+2+num_atom]]
                    return a, np.array([x,y,z]), c
                except (ValueError, IndexError):
                    raise ValueError(
                        f"Error parsing coordinate between lines {index_atom_start+2} and {index_atom_start+3+num_atom}"
                    )
        raise ValueError(
            "No ATOM section found"
        )

    @staticmethod
    def get_bonds(content: str|list) -> tuple:

        lines = read_lines(content, MOL2.file_extension) 
        num_bonds = MOL2.get_num_info(lines)[1]
        bonds, double_bonds, triple_bonds = [], [], []   

        for i, line in enumerate(lines): 
            if line.startswith("@<TRIPOS>BOND"):
                for line_bond in lines[i + 1 : i + 2 + num_bonds]:
                    tokens = line_bond.split()
                    a1, a2, bond_type = tokens[1], tokens[2], tokens[3]
                    # bond order is referenced as int values in str as from readlines 
                    if bond_type == "1":
                        bonds.append((a1, a2))
                    elif bond_type in ("2", "ar"):
                        double_bonds.append((a1, a2))
                    elif bond_type == "3":
                        triple_bonds.append((a1, a2))
                return bonds, double_bonds, triple_bonds
            else:
                raise ValueError(
                    "Cannot parse chemical information without MOLECULE, ATOM, BOND sections."
                )
        return bonds, double_bonds, triple_bonds

    @staticmethod
    def parse_mol2(content: str|list) -> dict:

        lines = read_lines(content, MOL2.file_extension) if isinstance(content, str) else content
        
        return {
            "num_atom"     : MOL2.get_num_info(lines)[0],
            "num_bond"     : MOL2.get_num_info(lines)[1],

            "mol_type"     : MOL2.get_mol_type(lines)[2],

            "charges_type" : MOL2.get_charges_types(lines)[2],

            "comment"      : MOL2.get_comment(lines),

            "atoms"        : MOL2.get_coord(lines)[0],
            "coords"       : MOL2.get_coord(lines)[1],
            "charges"      : MOL2.get_coord(lines)[2],

            "bonds"        : MOL2.get_bonds(lines)[0],
            "double_bonds" : MOL2.get_bonds(lines)[1],
            "triple_bonds" : MOL2.get_bonds(lines)[2],
        }

# %%
def parse_json(path: str):
    """
    JSON file do not have systematic representation.
    So they are just parsed using built in python parser.
    """

    with open(path, "r") as file:
        data = json.load(file)

    return data

# %%
def parse_yaml(path: str):
    """
    YAML file do not have systematic representation.
    So they are just parsed using built in function.
    """

    with open(path, "r") as file:
        data = yaml.safe_load_all(file)

    return data

# %%
class SDF:
    """
    Input supported :
        - Absolute path to file
        - List of read().strip() outputs
    """
    
    file_extension = "sdf"

    @staticmethod
    def read_sdf(content: str|list, extension: str) -> list:
        if isinstance(content, str):
            if content.split(".")[-1] == extension:
                with open(content, "r") as sdf:
                    list_molecules = sdf.read().strip()
                    return list_molecules.split("$$$$")
            else:
                raise ValueError(f"File is not {extension}")
        elif isinstance(content, list):
            return content
        else:
            raise TypeError("Input invalid Type")

    @staticmethod
    def get_num_info(content: str|list) -> list:
        blocks = SDF.parse_sdf(content)
        numinfos = []
        for block in blocks:
            lines = [l for l in block.strip().splitlines() if l.strip()]
            if lines is not None:
                num_info = MOL.get_num_info(lines)
                numinfos.append(num_info)
            else:
                numinfos.append([])
        return numinfos

    @staticmethod
    def get_properties(content: str|list) -> list:
        blocks = SDF.parse_sdf(content)
        mol_properties = []
        for block in blocks:
            lines = [l for l in block.strip().splitlines() if l.strip()]
            if lines is not None:
                properties = MOL.get_information(lines)
                mol_properties.append(properties)
            else:
                mol_properties.append([])
        return mol_properties
    
    @staticmethod
    def get_coord(content: str|list) -> list:
        blocks = SDF.parse_sdf(content)
        mol_coords = []
        for block in blocks:
            lines = [l for l in block.strip().splitlines() if l.strip()]
            if lines is not None:
                coords = MOL.get_coord(lines)
                mol_coords.append(coords)
            else:
                mol_coords.append([])
        return mol_coords

    @staticmethod    
    def get_bonds(content: str|list) -> list:
        blocks = SDF.parse_sdf(content)
        mol_bonds = []
        for block in blocks:
            lines = [l for l in block.strip().splitlines() if l.strip()]
            if lines is not None:
                bonds = MOL.get_bonds(lines)
                mol_bonds.append(bonds)
            else:
                mol_bonds.append([])
        return mol_bonds


    @staticmethod 
    def parse_sdf(content: str|list) -> list:
        """
        Parse SDF file into list of molecule dictionaries.
        Each molecule contains atoms, bonds, and property data.
        """

        blocks = SDF.read_sdf(content, SDF.file_extension)

        molecules = []
        for block in blocks:
            lines = [l for l in block.strip().splitlines() if l.strip()]
            if lines is not None:
                molecule = MOL.parse_mol(lines)
                molecules.append(molecule)
            else:                             
                molecules.append([])
        return molecules

# %%
class CML:

    @staticmethod
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

# %%
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

# %%
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
class Gaussian:

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
