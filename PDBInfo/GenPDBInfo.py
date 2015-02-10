'''
    Generate PDB information for each biounit file
'''
def runDSSP(model, pdbname):
    p = PDBParser()
    dssp = DSSP(model, pdbname)
    for residue in dssp:
        resinfo = residue[0]
        second_str = residue[1]
        ssa     = residue[2]
        rsa     = residue[3]


def EachBioUnit(biounit)
