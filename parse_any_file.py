# %%
import numpy as np
import json
import yaml

# %%
with open("/home/greg/Code/parse_any_chemical_file/ex.mol", "r") as mol:
    lines = mol.readlines()
pre_mol = list()
for i in range(1, int(lines[3].split()[0])):
    pre_mol.append(lines[3 + i].split()[0:4])
pre_mol = np.array(pre_mol)[:, [3, 0, 1, 2]]
print(pre_mol)


# %%
def check_len(molecule: dict) -> dict:
    """
    Check if molecule file is correct : testing if each coor col is correct : must be a np.ndarray
    of the following form :

    np.ndarray = | atom | x | y | z |
    If the order is not respected then the function cannot verrify that the given molecule is valid.

    Test made :

        - test if atoms are given with their symbol and not with their atomic number : should be modified with
        an existing periodic table
        - test if coords are only number is

    molecule (dict) : dictionnary obtained after using functions below
    """

    molecule_array = molecule["coord"]

    # testing atoms symbols
    test_array = np.vectorize(lambda x: isinstance(x, str))(molecule_array[:, 0])
    missing_pos = np.array([i for i in test_array.shape[0] if not test_array[:, 0]])

    if not any(test_array):
        raise ValueError(
            f"invalid file : atom at position {missing_pos + 1} is/are incorrect"
        )

    # testing coord numbers
    for col in range(1, molecule_array.shape[1]):
        test_array = np.vectorize(
            lambda x: isinstance(x, (int, float, complex, np.number))
        )(molecule_array[:, col])
        missing_pos = np.array([i for i in test_array.shape[0] if not test_array[i]])

        if not any(test_array):
            if col == 1:
                raise ValueError(
                    f"invalid file : coord x at position {missing_pos + 1} are incorrect"
                )

            elif col == 2:
                raise ValueError(
                    f"invalid file : coord x at position {missing_pos + 1} are incorrect"
                )

            elif col == 3:
                raise ValueError(
                    f"invalid file : coord x at position {missing_pos + 1} are incorrect"
                )

    return molecule


# %%
def parse_xyz(path: str) -> dict:
    """
    Open .xyz file and extract data as a dict
    Store molecule's coords as a np.ndarray in a single entry
    First and Second line are an entry.
    _____________________________________

    path (str) : path to the file to treat
    molecule (dict) : dict containing this structure :
    molecule = {
    "num_atom"
    "comment"
    "coord"
    }
    """

    # Read the file
    with open(path, "r") as xyz:
        lines = xyz.readlines()

    coords = list()
    num_atom = int(lines[0].split()[0])
    for line_coord in lines[2 : num_atom + 1]:
        tokens = line_coord.split()
        atoms = tokens[1]
        x, y, z = map(float, tokens[2:5])
        coords.append((atoms, x, y, z))

    # final object
    molecule = {
        "num_atom": num_atom,
        "comment": lines[1],
        "coord": np.array(coords),
    }

    return molecule


# %%
def parse_mol(path: str) -> dict:
    """
    Parse MOL like XYZ adding bonds and
    properties since its available.
    """

    # Read the file
    with open(path, "r") as mol:
        lines = mol.readlines()

    if "V3000" in lines[3].split():
        raise KeyError("V3000 formatting is not supported yet")

    # Get Coordinates information based on atoms number
    pre_mol = list()
    num_atom = int(lines[3].split()[0])
    for line_coord in lines[4 : num_atom + 1]:
        tokens = line_coord.split()
        atoms = tokens[3]
        x, y, z = map(float, tokens[0:3])
        pre_mol.append((atoms, x, y, z))

    # Get bonds
    bonds, double_bonds, triple_bonds = list(), list(), list()
    # TODO

    # Get optionnal information from the file
    j = 0
    properties = dict()
    line_info = lines[3 + j + num_atom]
    while not line_info.startswith("M END"):
        property = line_info.split()[1]
        values = line_info.split()[2:]
        properties[str(property)] = values
        j += 1

    molecule = {
        "num_atom": num_atom,
        "properties": properties,
        "bonds": bonds,
        "double_bonds": double_bonds,
        "triple_bonds": triple_bonds,
        "coord": np.array(pre_mol),
    }

    return molecule


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
            comment = lines[i + 5].strip() if len(lines) > i + 5 else None

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
    So they are just aprsed using built in function.
    """

    with open(path, "r") as file:
        data = yaml.safe_load_all(file)

    return data


def parse_sdf(path: str) -> list:
    """
    Parsing SDF files is very similar to .MOL since it's a container for multiple molecules under MOL format.
    Instead a list of dict (molecules) will be returned.
    """

    with open(path, "r") as sdf:
        lines = sdf.readlines()

    molecules = list()
    for line in lines:
        while not line.startswith("$$$$"):

            if "V3000" in lines[3].split():
                raise KeyError("V3000 formatting is not supported yet")

            # Get Coordinates information based on atoms number
            pre_mol = list()
            num_atom = int(lines[3].split()[0])
            for line_coord in lines[4 : num_atom + 1]:
                tokens = line_coord.split()
                atoms = tokens[3]
                x, y, z = map(float, tokens[0:3])
                pre_mol.append((atoms, x, y, z))

            # Get bonds
            bonds, double_bonds, triple_bonds = list(), list(), list()
            # TODO

            # Get optionnal information from the file
            j = 0
            properties = dict()
            line_info = lines[3 + j + num_atom]
            while not line_info.startswith("M END"):
                property = line_info.split()[1]
                values = line_info.split()[2:]
                properties[str(property)] = values
                j += 1

            molecule = {
                "num_atom": num_atom,
                "properties": properties,
                "bonds": bonds,
                "double_bonds": double_bonds,
                "triple_bonds": triple_bonds,
                "coord": np.array(pre_mol),
            }
            molecules.append(molecule)

    return molecules


def parse_mdl(path: str) -> dict:
    """
    Parse MDL like
    """

    tree = ET.parse(path)
    root = tree.getroot()
    ns = {'cml': 'http://www.xml-cml.org/schema'}

    atoms = [atom.get("elementType") for atom in root.findall(".//cml:atom", ns)]
    if 
    x = [atom.get("x3") for atom in root.findall(".//cml:atom", ns)]
    y = [atom.get("y3") for atom in root.findall(".//cml:atom", ns)]
    z = [atom.get("z3") for atom in root.findall(".//cml:atom", ns)]

    return molecule
