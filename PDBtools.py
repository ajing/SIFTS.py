'''
    Tools to retrievel info from PDB
'''

import sys
import os
from Bio.PDB import PDBParser
import gzip
from SCOPData import protein_letters_3to1

#DISABLE_RM_IN_DICT = True
DISABLE_RM_IN_DICT = False
#BIOLIP_DIR = True
BIOLIP_DIR = False

if BIOLIP_DIR:
    from AABindingSiteDist.Model import PDBDIR
else:
    from SIFTS.SIFTSXMLMapModel import PDBDIR

def GetFilewithPDB(pdbid):
    if BIOLIP_DIR:
        pdbfname = pdbid + ".pdb.gz"
        return os.path.join(PDBDIR, pdbfname)

    pdbfname = "pdb" + pdbid + ".ent.gz"
    subdir   = pdbid[1:3]
    return os.path.join(PDBDIR, subdir, pdbfname)

def PrintPDBstruct(pdbstruct):
    for model in pdbstruct:
        for chain in model:
            for residue in chain:
                print residue

def GetResidueFromPDB(pdbstruct, chainid = None, resname = None, resnum = None):
    if chainid is None and resname is None and resnum is None:
        reslist = []
        for res in pdbstruct[0].get_residues():
            if res.get_resname() in protein_letters_3to1:
                reslist.append(res)
        return reslist
    elif resname is None and resnum is None:
        reslist = []
        for res in pdbstruct[0][chainid]:
            if res.get_resname() in protein_letters_3to1:
                reslist.append(res)
        return reslist

    for model in pdbstruct:
        try:
            residue = model[chainid][int(resnum)]
            ### for simplicity
            return residue
            #if residue.get_resname() == resname:
            #    return residue
        except:
            continue
    raise Exception("Cannot find chain: %s, resname: %s, resnum: %s to PDB:%s" % ( chainid, resname, resnum, pdbstruct))

def GetStructure(pdbid):
    parser = PDBParser(PERMISSIVE = 1)
    pdbobj = gzip.open(GetFilewithPDB(pdbid))
    structure = parser.get_structure(pdbid, pdbobj)
    return structure

pdbdict = dict()
def GetResidueObj(pdbid, chainid = None, resname = None, resnum = None):
    if pdbid in pdbdict:
        struct = pdbdict[pdbid]
    else:
        struct  = GetStructure(pdbid)
        if not DISABLE_RM_IN_DICT:
            for eachkey in pdbdict.keys():
                if eachkey != pdbid:
                    del pdbdict[eachkey]
        pdbdict[pdbid] = struct
    if chainid is None and resname is None and resnum is None:
        return GetResidueFromPDB(struct)
    elif resname is None and resnum is None:
        return GetResidueFromPDB(struct, chainid)
    return GetResidueFromPDB(struct, chainid, resname, resnum)


# copy and gunzip file
def CopyAndGunzip(infiledir, outfiledir):
    inF = gzip.open(infiledir, 'rb')
    outF = open(outfiledir, 'wb')
    outF.write( inF.read() )
    inF.close()
    outF.close()

if __name__ == "__main__":
    reslist = GetResidueObj("2ml1")
    for res in reslist:
        print res
