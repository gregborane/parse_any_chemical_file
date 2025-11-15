# %%
import sys
import json
from itertools import combinations
from scipy.spatial.distance import euclidean
import numpy as np

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

    try : 
        atom1_3dcoord = coord[atom1,:]
        atom2_3dcoord = coord[atom2,:]
    except Exception: 
        raise ValueError("Invalid indices given to function")

    return np.float64(euclidean(atom1_3dcoord, atom2_3dcoord))


# %%
def generate_bonds(content: str, tolerance: np.float16 = np.float16(0.1)):
    global periodic_table, list_p_table_keys

    coord, atoms, bonds, double_bonds, triple_bonds = get_molecule(content)
    def find_bonds(bonds_type : str) -> set:
        
        existing_bonds = set()
        for i, j in combinations(range(len(atoms)), 2):
            elem1 = atoms[i] if isinstance(atoms[i], str) else list_p_table_keys[atoms[i]]
            elem2 = atoms[j] if isinstance(atoms[j], str) else list_p_table_keys[atoms[j]]
            dist = compute_distance(i, j, coord)
            
            # l(bond) ~ r(Elem1) + r(Elem2)
            radius_sum = periodic_table[elem1][bonds_type] + periodic_table[elem2][bonds_type]
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

    

    return bonds, double_bonds, triple_bonds

