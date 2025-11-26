# %%
import sys
import json
from itertools import combinations
from scipy.spatial.distance import euclidean
import numpy as np
from collections import defaultdict
from rkdit import Chem
import networkx as nx

sys.path.append("./")

with open("./constant/covalences.json", "r") as p_table:
    periodic_table = json.load(p_table)
list_p_table_keys = list(periodic_table.keys())

print(periodic_table[list_p_table_keys[0]])


# %%
def get_molecule(content: str | list) -> tuple:
    from parse_any_file import parse_any_file

    coord = parse_any_file(content)["coord"]
    atoms = parse_any_file(content)["atoms"]
    bonds = parse_any_file(content)["bonds"]
    double_bonds = parse_any_file(content)["double_bonds"]
    triple_bonds = parse_any_file(content)["triple_bonds"]

    return coord, atoms, bonds, double_bonds, triple_bonds


# %%
def compute_distance(
    atom1: int,
    atom2: int,
    coord: np.ndarray,
) -> np.float64:
    try:
        atom1_3dcoord = coord[atom1, :]
        atom2_3dcoord = coord[atom2, :]
    except Exception:
        raise ValueError("Invalid indices given to function")

    return np.float64(euclidean(atom1_3dcoord, atom2_3dcoord))


# %%
def generate_bonds(content: str, tolerance: np.float16 = np.float16(0.1)):
    global periodic_table, list_p_table_keys

    coord, atoms, bonds, double_bonds, triple_bonds = get_molecule(content)

    def find_bonds(bonds_type: str) -> set:
        existing_bonds = set()
        for i, j in combinations(range(len(atoms)), 2):
            elem1 = (
                atoms[i] if isinstance(atoms[i], str) else list_p_table_keys[atoms[i]]
            )
            elem2 = (
                atoms[j] if isinstance(atoms[j], str) else list_p_table_keys[atoms[j]]
            )
            dist = compute_distance(i, j, coord)

            # l(bond) ~ r(Elem1) + r(Elem2)
            radius_sum = (
                periodic_table[elem1][bonds_type] + periodic_table[elem2][bonds_type]
            )
            if dist < radius_sum + tolerance:
                existing_bonds.add((i, j))
        return existing_bonds

    if not bonds:
        bonds = find_bonds("single_covalent_radii")
    if not double_bonds:
        double_bonds = find_bonds("double_covalent_radii")
    if not triple_bonds:
        triple_bonds = find_bonds("triple_covalent_radii")

    """
    Conjugated system will give extra bonds because bonds order is not an int which is the case in our consideration.

    To correct such, we need to algorithmically remove those extra bonds by checking valencies of all present atoms. 
    H = 1, C = 4, etc ... 

    Such determination depends on the valencies considered for the atoms. Hypervalencies atoms will not be taken into account. 
    """

    # Step 1: Count current bonds for each atom
    bond_count = defaultdict(int)
    for i, j in bonds:
        bond_count[i] += 1
        bond_count[j] += 1

    # Step 1: Remove unecessayry bonds from mono valent atoms
    for i, j in bonds.copy():
        if (atoms[i] in ["H", "F", "Cl", "Br", "I"] and atoms[j] in ["W", "Mo"]) or (
            atoms[j] in ["H", "F", "Cl", "Br", "I"] and atoms[i] in ["W", "Mo"]
        ):
            bonds.remove((i, j))
            bond_count[i] -= 1
            bond_count[j] -= 1

    # Set 2: Get Mo=C bonds
    if not check_cond:
        for i, j in bonds:
            if atoms[i] == "C" and atoms[j] in ["W", "Mo"]:
                if (
                    bond_count[i] < valencies["C"]
                    and bond_count[j] < valencies[atoms[j]]
                ):  # Ensure valency allows it
                    if (i, j) in bonds:
                        double_bonds.add((i, j))  # Convert single bond to double
                        bond_count[i] += 1
                        bond_count[j] += 1

            elif atoms[j] == "C" and atoms[i] in ["W", "Mo"]:
                if (
                    bond_count[j] < valencies["C"]
                    and bond_count[i] < valencies[atoms[i]]
                ):  # Ensure valency allows it
                    if (i, j) in bonds:
                        double_bonds.add((i, j))  # Convert single bond to double
                        bond_count[i] += 1
                        bond_count[j] += 1
    else:
        double_bonds.add((1, 22))  # Convert single bond to double
        bonds.add((32, 29))
        bonds.remove((22, 29))

        bond_count[22] += 1
        bond_count[1] += 1
        bond_count[32] += 1

    # Step 3: Ensure Carbon (C) has exactly 4 bonds
    for i, atom in enumerate(atoms):
        if atom == "C" and bond_count[i] == 3:
            for j in list(b for a, b in bonds if a == i) + list(
                a for a, b in bonds if b == i
            ):
                if (
                    bond_count[j] < valencies.get(atoms[j], 4)
                    and (i, j) not in double_bonds
                ):
                    double_bonds.add((i, j))  # Convert single bond to double
                    bond_count[i] += 1
                    bond_count[j] += 1
                    break

    return bonds, double_bonds, triple_bonds


def get_graph(cat: dict) -> nx.Graph:
    """
    Convert xyz file to a nx.graph object
    bonds are edges, nodes are atoms
    cat (dict) : dict containing a four coordinates object with atoms and its respective coordinates
    return G (nx.Graph) : graph object of the molecule
    """
    G = nx.Graph()
    # Add nodes to the graph based on atom positions and symbols
    for i, atom in enumerate(cat["atoms"]):
        G.add_node(i, element=atom, pos=(cat["x"][i], cat["y"][i], cat["z"][i]))

    # Add edges for single bonds
    for bond in cat["bonds"]:
        if bond in cat["double_bonds"]:
            u, v = bond
            bond_type = (
                "double"  # Testing if this bonds has been identified as double already
            )
            G.add_edge(u, v, bond_type=bond_type)

        else:
            u, v = bond
            bond_type = "single"  # Attributing the single type is no other exist
            G.add_edge(u, v, bond_type=bond_type)
    return G


def nx_to_mol(G: nx.Graph):
    """
    Convert a networkx graph boject to mol one
    G (nx.graph): graph to convert
    return mol : rdkit mol object
    """
    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, "element")
    # atomic_nums = [periodic_table[atom] for atom in atoms.values()]
    node_to_idx = {}
    for node in G.nodes():
        a = Chem.Atom(atomic_nums[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    bond_types = nx.get_edge_attributes(G, "bond_type")
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type_map[bond_type])

    Chem.SanitizeMol(mol)
    return mol
