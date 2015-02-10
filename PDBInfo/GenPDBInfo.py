'''
    Generate PDB information for each biounit file
'''

import os
from Bio.PDB import PDBParser
from Bio.PDB import DSSP

OUTOBJ = "dsspout_pdb.txt"
BIODIR = "biounit"

def RunDSSP(model, pdbfile):
    dssp = DSSP(model, pdbfile)
    reslist = []
    for residue in dssp:
        resinfo = residue[0]
        second_str = residue[1]
        ssa     = residue[2]
        rsa     = residue[3]
        reslist.append({"res_obj": resinfo, "sec_str": second_str, "ssa": ssa, "rsa": rsa})
    return reslist

def ProcessDSSP(reslist):
    newlist = []
    for eachres in reslist:
        residue = eachres["res_obj"]
        resid = residue.get_full_id()
        # PDBID, model id, chain id, residue name, residue num, secondary structure, ssa, rsa
        newlist.append([resid[0], resid[1], resid[2], residue.resname, resid[3][1], eachres["sec_str"], eachres["ssa"], eachres["rsa"]])
    return newlist

def RunEachBioUnit(biounit):
    p = PDBParser()
    pdbname= biounit.split("/")[-1]
    models = p.get_structure(pdbname, biounit)
    outlines = []
    for model in models:
        dssp_model = RunDSSP(model, biounit)
        lines      = ProcessDSSP(dssp_model)
        outlines  += lines
    return outlines

def AllBioUnit(directory):
    out_obj = open(OUTOBJ, "w")
    for eachbiounit in os.listdir(directory):
        biounit = os.path.join(directory, eachbiounit)
        lines   = RunEachBioUnit(biounit)
        if lines:
            content = "\n".join(["\t".join(map(str, x)) for x in lines])
            out_obj.write(content)

if __name__ == "__main__":
    #EachBioUnit("pdb10gs.ent")
    AllBioUnit(BIODIR)
