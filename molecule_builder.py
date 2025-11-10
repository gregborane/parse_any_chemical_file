import sys, os
import networkx

sys.path.append("./")

from parse_any_file import parse_any_file


class BuildRdkit:
    def __init__(self, content: str | list):
        self.coord = parse_any_file(content)["coord"]
        self.atoms = parse_any_file(content)["atoms"]
        self.bonds = parse_any_file(content)["bonds"]
        self.double_bonds = parse_any_file(content)["double_bonds"]
        self.triple_bonds = parse_any_file(content)["triple_bonds"]

    # Magic number =
    ## Initial object

    def generate_bonds(self):
        if all((self.bonds, self.double_bonds, self.triple_bonds)):
            return self.bonds, self.double_bonds, self.triple_bonds

        return self.bonds, self.double_bonds, self.triple_bonds

    def create_graph(self):
        return self.atoms, self.coord
