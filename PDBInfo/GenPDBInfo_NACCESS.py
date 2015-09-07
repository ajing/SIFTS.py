'''
    Generate PDB information for each biounit file
'''

import os
from Bio.PDB import PDBParser
from NACCESS import NACCESS

BIODIR = "2013_biounits_noligand"
#BIODIR = "/home/ajing/Documents/Research/SIFTS.py/PDBInfo/test"
BIODIR = "../../ligandNet/2013_biounits_noligand"
OUTDIR = "/home/ajing/Documents/Research/SIFTS.py/PDBInfo/out_naccess"
NACCESS_DIR = "/home/ajing/Documents/Research/SIFTS.py/PDBInfo/naccess2.1.1/naccess"
TMP_DIR = "/home/ajing/Documents/Research/SIFTS.py/PDBInfo/tmp"

COLNAME = ["all_atoms_abs", "all_atoms_rel", "non_polar_abs", "non_polar_rel", "all_polar_abs", "all_polar_rel"]

def RunNACCESS(model, pdbfile):
    try:
        naccess = NACCESS(model, pdbfile, NACCESS_DIR, TMP_DIR)
    except:
        return None
    reslist = []
    for residue in naccess:
        resinfo = residue[0]
        reslist.append(dict({"res_obj": resinfo}.items() + residue[1].items()))
    return reslist

def ProcessDSSP(reslist):
    newlist = []
    for eachres in reslist:
        residue = eachres["res_obj"]
        resid = residue.get_full_id()
        # PDBID, model id, chain id, residue name, residue num, secondary structure, ssa, rsa
        newline = [resid[0].split(".")[-2][-4:], resid[0], resid[1], resid[2], residue.resname, resid[3][1]]
        for names in COLNAME:
            newline.append(eachres[names])
        newlist.append(newline)
    return newlist

def RunEachBioUnit(biounit):
    p = PDBParser(PERMISSIVE = 1)
    pdbname= biounit.split("/")[-1]
    try:
        #print "models for:", pdbname, biounit
        models = p.get_structure(pdbname, biounit)
    except:
        return None
    outlines = []
    for model in models:
        dssp_model = RunNACCESS(model, biounit)
        if dssp_model:
            lines      = ProcessDSSP(dssp_model)
            outlines  += lines
    return outlines

#def FileFilter(filelist, exist_dir):
#    outfiles = os.listdir(exist_dir)
#    outfileb = [x.split(".")[:-4] for x in outfiles]
#    return [x for x in filelist if not x in outfileb]

def FileFilter(filelist, exist_dir):
    return filelist

def AllBioUnit(directory):
    fileleft = FileFilter(os.listdir(directory), OUTDIR)
    for eachbiounit in fileleft:
        biounit = os.path.join(directory, eachbiounit)
        output  = eachbiounit + ".out"
        outobj  = open(os.path.join(OUTDIR, output), "w")
        lines   = RunEachBioUnit(biounit)
        if lines:
            content = "\n".join(["\t".join(map(str, x)) for x in lines])
            outobj.write(content)
        outobj.close()

if __name__ == "__main__":
   # RunEachBioUnit("10gs.bio1")
    AllBioUnit(BIODIR)
