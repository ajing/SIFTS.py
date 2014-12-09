'''
    Get the residue depth for each residue in BioLiP
    run as:
    python -m ResidueDepth.Controller
'''
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.ResidueDepth import get_surface, residue_depth, ca_depth
from Bio.PDB.Polypeptide import is_aa
import os
from AABindingSiteDist.Controller import BSParser
from PDBtools import GetFilewithPDB, CopyAndGunzip, GetStructure

PDBTOXYZ = "./ResidueDepth/msms/pdb_to_xyzr"
MSMS = "./ResidueDepth/msms/msms.x86_64Linux2.2.6.1"
OUTCA = "aveResCaDep.txt"
OUTALL = "aveResAllDep.txt"
OUT   = "avedist2surface.txt"
BIOLIP_DIR = "./Data/bindingsite2.txt"

# working directory
WDIR  = "./ResidueDepth/tmp"

def GetResidueDepPDB(pdb, pdbfile):
    s  = GetStructure(pdb)
    model = s[0]
    residuelist = Selection.unfold_entities(model, 'R')
    try:
        surface = get_surface(pdbfile, PDBTOXYZ, MSMS)
    except:
        print "cannot get surface for " + pdbfile
        return
    outobj = open(OUT, "a")
    for residue in residuelist:
        if not is_aa(residue):
            continue
        # minimun average depth for all atoms
        resid   = residue.get_id()
        resname = residue.get_resname()
        chainid = residue.get_parent().get_id()
        try:
            rd = residue_depth(residue, surface)
        except:
            continue
        ca_rd = ca_depth(residue, surface)
        info    = [pdb, chainid, resid[1], resname, str(rd), str(ca_rd)]
        #print info
        outobj.write("\t".join(map(str, info)) + "\n")
    outobj.close()

def RunAllBioLiPPDB():
    bslist = BSParser(BIOLIP_DIR)
    pdblist = []
    try:
        os.remove(OUT)
    except:
        pass
    for bs in bslist:
        pdb = bs.pdbid
        if not pdb in pdblist:
            pdblist.append(pdb)
    for pdb in pdblist:
        outdir = os.path.join(WDIR, pdb)
        pdbfile = GetFilewithPDB(pdb)
        CopyAndGunzip(pdbfile, outdir)
        GetResidueDepPDB(pdb, outdir)

if __name__ == "__main__":
    pdbfile = "./tmp/pdb110m.ent"
    #GetResidueDepPDB("110m", pdbfile)
    RunAllBioLiPPDB()
