'''
   testing for NACCESS
'''

from Bio.PDB import PDBParser
from NACCESS import NACCESS

NACCESS_DIR = "/home/ajing/Documents/Research/SIFTS.py/PDBInfo/naccess2.1.1/naccess"
TMP_DIR     = "tmp"

def TestNACCESS(model, pdbfile):
    naccess = NACCESS(model, pdbfile, NACCESS_DIR)
    for residue in naccess:
        print residue

def main():
    p = PDBParser()
    filename = "test/10gs.bio1"
    models = p.get_structure("10gs", filename)
    for model in models:
        print models[0]
        print model.get_full_id()
        TestNACCESS(models[0], filename)

if __name__ == "__main__":
    main()
