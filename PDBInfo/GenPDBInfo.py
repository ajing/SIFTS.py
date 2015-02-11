'''
    Generate PDB information for each biounit file
'''

import os
from Bio.PDB import PDBParser
from Bio.PDB import DSSP

OUTOBJ = "dsspout_pdb.txt"
BIODIR = "2013_biounits"
OUTDIR = "out"

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
        newlist.append([resid[0].split(".")[-2][-4:], resid[0], resid[1], resid[2], residue.resname, resid[3][1], eachres["sec_str"], eachres["ssa"], eachres["rsa"]])
    return newlist

def RunEachBioUnit(biounit):
    p = PDBParser(PERMISSIVE = 1)
    pdbname= biounit.split("/")[-1]
    try:
        models = p.get_structure(pdbname, biounit)
    except:
        return None
    outlines = []
    for model in models:
        dssp_model = RunDSSP(model, biounit)
        lines      = ProcessDSSP(dssp_model)
        outlines  += lines
    return outlines

def FileFilter(filelist):
    outfiles = os.listdir(OUTDIR)
    outfileb = [x.split(".")[-2] for x in outfiles]
    return [x for x in filelist if not x in outfileb]

def AllBioUnit(directory):
    for eachbiounit in FileFilter(os.listdir(directory)):
        biounit = os.path.join(directory, eachbiounit)
        output  = eachbiounit + ".out"
        outobj  = open(os.path.join(OUTDIR, output), "w")
        lines   = RunEachBioUnit(biounit)
        if lines:
            content = "\n".join(["\t".join(map(str, x)) for x in lines])
            outobj.write(content)
        outobj.close()

if __name__ == "__main__":
    #EachBioUnit("pdb10gs.ent")
    AllBioUnit(BIODIR)
