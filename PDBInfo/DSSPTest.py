'''
   Just for testing
'''

from Bio.PDB import PDBParser
from Bio.PDB import DSSP

def TestDSSP(model, pdbfile):
    dssp = DSSP(model, pdbfile)
    resinfo = []
    #for res in model.get_residues():
    #    print res
    for residue in dssp:
        resinfo = residue[0]
        second_str = residue[1]
        ssa     = residue[2]
        rsa     = residue[3]
    print dir(resinfo)
    print resinfo.get_full_id()
    print resinfo.)

def main():
    p = PDBParser()
    filename = "pdb10gs.ent"
    models = p.get_structure("10GS", filename)
    for model in models:
        print models[0]
        print model.get_full_id()
        TestDSSP(models[0], filename)

if __name__ == "__main__":
    main()
