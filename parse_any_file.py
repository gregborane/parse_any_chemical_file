# %%
import numpy as np

# %%
def parse_xyz(path):

    with open(path, "r") as xyz:
        lines = xyz.readlines()

    atom_coord = np.array(lines[2:])
    actual_len = len(atom_coord)

    i = 0
    for col in range(atom_coord.shape[1]):

        if len(~np.isnan(atom_coord[:,col])) != actual_len:
            return "invalid file not all lines have the same size"

        else:
            i += 1

    if i == atom_coord.shape[0]:
        print("valid file")
        return atom_coord


# %%
