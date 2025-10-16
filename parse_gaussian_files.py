# %%
import numpy as np
import json
import yaml
import xml.etree.ElementTree as ET
from rdkit import Chem

class ReadGaussian:
    def __init__(self, path: str):
       self.path = path

    def test_extensions(self):
        """
        Verify validity of the extension given
        """
        path_file = self.path

        if path_file.split(".")[-1] not in ["cube", "gjf", "log"]:
            raise ValueError("Extensions is not a valid gaussian file")
    
        return True

    def read_lines(self):
        """
        Get Data
        """
        if self.test_extensions(): 
            with open(self.path,"r") as gaussian_file:
                    lines = gaussian_file.readlines()
                    self.lines = lines
        return self.lines

    def extract_coordinates(self):
        """
        Obtain Atoms coordinates and atomic number
        """
        lines = self.read_lines()
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

    def extract_charges(self, path: str):
        """
        Obtain partial charges from NBO calculations
        Mulliken, 
        """

        lines = self.read_lines()
        indice_mulliken = None, 

        for i, lines in enumerate(lines[::-1]):
            split_line = line.split()

            if "Mulliken" in lines.split() and len(lines.split()) == 2:
                indice_mulliken = i

            if "APT" in lines.split() and len(lines.split()) == 2:

        return "bite"

    def extract_bond_order(self, path: str):
        """
        Obtain bond order from NBO calculation
        """
        return 'bite bite'

# %%
listm = "NAtoms=    9 NActive=    9 NUniq=    9 SFac= 1.00D+00 NAtFMM=   60 NAOKFM=F Big=F"
print(listm.split())

# %%

