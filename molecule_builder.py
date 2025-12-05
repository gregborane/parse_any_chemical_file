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

# %%
with open("./constant/covalences.json", "r") as p_table:
    periodic_table = json.load(p_table)
list_p_table_keys = list(periodic_table.keys())


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
def generate_bonds(
    content: str | list, tolerance: np.float16 = np.float16(0.1)
) -> tuple:
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
    Most sensitive atom is Carbon because it is often implied in aromatic bonding 

    Such determination depends on the valencies considered for the atoms. Hypervalencies atoms will not be taken into account. 
    """

    bond_count = defaultdict(int)
    for i, j in bonds:
        bond_count[i] += 1
        bond_count[j] += 1

    for i in range(len(atoms)):
        elem1 = list_p_table_keys[atoms[i]]
        if elem1 == "C" and bond_count[i] == 3:
            for j in list(b for a, b in bonds if a == i) + list(
                a for a, b in bonds if b == i
            ):
                valencies = periodic_table[list_p_table_keys[j]]["covalence"]
                if bond_count[j] < valencies and (i, j) not in double_bonds:
                    double_bonds.add((i, j))  # Convert single bond to double
                    bond_count[i] += 1
                    bond_count[j] += 1
                    break

    return bonds, double_bonds, triple_bonds


def get_graph(content: str | list) -> nx.Graph:
    coord, atoms, bonds, double_bonds, triple_bonds = get_molecule(content)
    bonds, double_bonds, triple_bonds = generate_bonds(content)
    G = nx.Graph()
    # Add nodes to the graph based on atom positions and symbols
    for i in range(len(atoms)):
        elem = periodic_table[list_p_table_keys[i]]
        G.add_node(i, element=elem, pos=(coord[0][i], coord[1][i], coord[2][i]))

    # Add edges for single bonds
    for bond in bonds:
        if bond in triple_bonds:
            u, v = bond
            bond_type = "triple"
            G.add_edge(u, v, bond_type=bond_type)
        elif bond in double_bonds:
            u, v = bond
            bond_type = "double"
            G.add_edge(u, v, bond_type=bond_type)
        else:
            u, v = bond
            bond_type = "single"
            G.add_edge(u, v, bond_type=bond_type)
    return G


def nx_to_mol(G: nx.Graph) -> Chem.RWMol:
    bond_type_map = {
        "single": Chem.BondType.SINGLE,
        "double": Chem.BondType.DOUBLE,
        "triple": Chem.BondType.TRIPLE,
    }

    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, "element")
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


def main(content: str | list) -> dict:
    coord, atoms, bonds, double_bonds, triple_bonds = get_molecule(content)
    bonds, double_bonds, triple_bonds = generate_bonds(content)

    graph_mol = get_graph(content)
    rdkit_mol = nx_to_mol(graph_mol)

    return {
        "coord": coord,
        "atoms": atoms,
        "bonds": bonds,
        "double_bonds": double_bonds,
        "triple_bonds": triple_bonds,
        "graph_mol": graph_mol,
        "rdkit_mol": rdkit_mol,
    }
