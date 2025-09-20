# %%
import numpy as np

# %%
def check_len(molecule_array : np.ndarray) -> np.ndarray:

    """
    Check if molecule file is correct : testing if each col is correct.
        - test if atoms are given with their symbol and not with their atomic number 
        - test if coords are only number is 
    molecule_array (np.array) : array of size 4*num_atom : atom | x | y | z | 
    num_atom : integer representing the number of atom in the molecule
    """
    
    # testing atoms symbols
    test_array = np.vectorize(lambda x: isinstance(x, str))(molecule_array[:,0])
    missing_pos = np.array([i for i in test_array.shape[0] if not test_array[:,0]]) 
    
    if not any(test_array):
        raise ValueError(f"invalid file : atom at position {missing_pos + 1} is/are incorrect") 

    # testing coord numbers
    for col in range(1, molecule_array.shape[1]):
        test_array = np.vectorize(lambda x: isinstance(x, (int, float, complex, np.number)))(molecule_array[:,col])
        missing_pos = np.array([i for i in test_array.shape[0] if not test_array[i]]) 
        
        if not any(test_array) or molecule_array:
            if col == 1:
                raise ValueError(f"invalid file : coord x at position {missing_pos + 1} are incorrect")
            
            elif col == 2:
                raise ValueError(f"invalid file : coord x at position {missing_pos + 1} are incorrect")
            
            elif col == 3:
                raise ValueError(f"invalid file : coord x at position {missing_pos + 1} are incorrect")
            
    return molecule_array

# %%
def parse_xyz(path : str) -> np.ndarray :
    """
    Open .xyz file and extract the four columns counting from the 3rd line 
    where start to appear info are
    """

    with open(path, "r") as xyz:
        lines = xyz.readlines()

    return np.array(lines[2:])

# %%
def parse_mol(path : str) -> np.ndarray :
    return 12


