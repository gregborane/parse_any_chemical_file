# %%
import numpy as np

# %%
with open("/home/greg/Code/parse_any_chemical_file/ex.mol", "r") as mol :
    lines = mol.readlines()
print(lines[3].split()[0])

# %%
def check_len(molecule : dict) -> dict :

    """
    Check if molecule file is correct : testing if each col is correct.
        - test if atoms are given with their symbol and not with their atomic number 
        - test if coords are only number is 
    molecule (dict) : array of size 4*num_atom : atom | x | y | z | 
    """

    molecule_array = molecule["coord"] 
    
    # testing atoms symbols
    test_array = np.vectorize(lambda x: isinstance(x, str))(molecule_array[:,0])
    missing_pos = np.array([i for i in test_array.shape[0] if not test_array[:,0]]) 
    
    if not any(test_array):
        raise ValueError(f"invalid file : atom at position {missing_pos + 1} is/are incorrect") 

    # testing coord numbers
    for col in range(1, molecule_array.shape[1]):
        test_array = np.vectorize(lambda x: isinstance(x, (int, float, complex, np.number)))(molecule_array[:,col])
        missing_pos = np.array([i for i in test_array.shape[0] if not test_array[i]]) 
        
        if not any(test_array):
            if col == 1:
                raise ValueError(f"invalid file : coord x at position {missing_pos + 1} are incorrect")
            
            elif col == 2:
                raise ValueError(f"invalid file : coord x at position {missing_pos + 1} are incorrect")
            
            elif col == 3:
                raise ValueError(f"invalid file : coord x at position {missing_pos + 1} are incorrect")
            
    return molecule

# %%
def parse_xyz(path : str) -> dict :
    """
    Open .xyz file and extract data as a dict
    Store molecule's coords as a np.ndarray in a single entry
    Rest is a dict entry
    """

    with open(path, "r") as xyz:
        lines = xyz.readlines()
    
    molecule = {
            "num_atom"   : lines[0], 
            "bonus info" : lines[1],
            "coord"      : np.array(lines[2:])
    }

    return molecule

# %%
def parse_mol(path : str) -> dict :
    """
    Open .mol file and extract data as dict
    Store molecule's coords as a np.ndarray
    Rest is a dict entry
    """

    with open(path, "r") as mol:
        lines = mol.readlines()

    pre_mol = list()
    for i in range(len(lines[3].split()[0])):
        pre_mol.append(lines[3+i].split()[0:3])
    
    pre_mol = np.array(pre_mol)[:,[3,0,1,2]]
    molecule = {
            "num_atom"   : lines[3].split()[0], 
            "bonus info" : "skipped_atm",
            "coord"      : pre_mol
    }

    return molecule


