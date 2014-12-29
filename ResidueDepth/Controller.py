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

from multiprocessing import Pool
import threading

DEBUG  = False

PDBTOXYZ = "./ResidueDepth/msms/pdb_to_xyzr"
MSMS = "./ResidueDepth/msms/msms.x86_64Linux2.2.6.1"
OUTCA = "aveResCaDep.txt"
OUTALL = "aveResAllDep.txt"
OUT   = "avedist2surface.txt"
BIOLIP_DIR = "./Data/bindingsite2.txt"

if DEBUG:
    OUTCA = OUTCA + "_tmp"
    OUTALL = OUTALL + "_tmp"
    OUT   = OUT + "_tmp"

# working directory
WDIR  = "./ResidueDepth/tmp"
# thread safe for writing file
mutex_writefile = threading.Lock()

def GetResidueDepPDB(pdb, pdbfile):
    s  = GetStructure(pdb)
    model = s[0]
    residuelist = Selection.unfold_entities(model, 'R')
    try:
        surface = get_surface(pdbfile, PDBTOXYZ, MSMS)
    except:
        print "cannot get surface for " + pdbfile
        return
    content = ""
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
        for each in info:
            if not each:
                continue
        #print info
        newline = "\t".join(map(str, info)) + "\n"
        content = content + newline

    mutex_writefile.acquire()
    outobj = open(OUT, "a")
    outobj.write(content)
    outobj.close()
    mutex_writefile.release()

def RemoveExistingPDB(pdblist):
    existpdbs  = []
    newpdblist = []
    for line in open(OUT):
        content = line.split()
        pdb     = content[0]
        if not pdb in existpdbs:
            existpdbs.append(pdb)
    print len(existpdbs)
    for eachpdb in pdblist:
        if not eachpdb in existpdbs:
            newpdblist.append(eachpdb)
    print len(newpdblist)
    return newpdblist

def RunOnePDB(pdb):
    outdir = os.path.join(WDIR, pdb)
    pdbfile = GetFilewithPDB(pdb)
    CopyAndGunzip(pdbfile, outdir)
    GetResidueDepPDB(pdb, outdir)

def RunAllBioLiPPDB():
    bslist = BSParser(BIOLIP_DIR)
    pdblist = []
    #try:
    #    os.remove(OUT)
    #except:
    #    pass
    for bs in bslist:
        pdb = bs.pdbid
        if not pdb in pdblist:
            pdblist.append(pdb)
    print "Number of PDBs before remove existing PDBs:", len(pdblist)
    pdblist = RemoveExistingPDB(pdblist)
    print "Number of PDBs after remove existing PDBs:", len(pdblist)
    print "one example:", pdblist[0]
    #for pdb in pdblist:
    #    print pdb

    pool = Pool(processes = 4)
    result = pool.map_async( RunOnePDB, pdblist)
    resulttxt = result.wait()
    print resulttxt

if __name__ == "__main__":
    pdbfile = "./tmp/pdb110m.ent"
    #GetResidueDepPDB("110m", pdbfile)
    #RemoveExistingPDB("")
    RunAllBioLiPPDB()
