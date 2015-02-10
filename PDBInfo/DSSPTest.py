'''
   Just for testing
'''

from Bio.PDB import PDBParser
from Bio.PDB import DSSP

def TestDSSP(model, pdbname):
    dssp = DSSP(model, pdbname)
    for residue in dssp:
        resinfo = residue[0]
        second_str = residue[1]
        ssa     = residue[2]
        rsa     = residue[3]

def main():
    p = PDBParser()
    model = p.get_structure("10GS", "pdb10gs.ent")[0]
    TestDSSP(model, "10GS")

if __name__ == "__main__":
    main()
