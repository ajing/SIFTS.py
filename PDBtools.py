'''
    Tools to retrievel info from PDB
'''

from SIFTSXMLMapModel import PDBDIR
import os
from Bio.PDB import PDBParser
import gzip

def GetFilewithPDB(pdbid):
    pdbfname = "pdb" + pdbid + ".ent.gz"
    return os.path.join(PDBDIR, pdbfname)

def PrintPDBstruct(pdbstruct):
    for model in pdbstruct:
        for chain in model:
            for residue in chain:
                print residue

def GetResidueFromPDB(pdbstruct, chainid, resname, resnum):
    for model in pdbstruct:
        try:
            residue = model[chainid][int(resnum)]
            if residue.get_resname() == resname:
                return residue
        except:
            continue
    raise Exception("Cannot find chain: %s, resname: %s, resnum: %s to PDB:%s" % ( chainid, resname, resnum, pdbstruct))

def GetStructure(pdbid):
    parser = PDBParser(PERMISSIVE = 1)
    pdbobj = gzip.open(GetFilewithPDB(pdbid))
    structure = parser.get_structure(pdbid, pdbobj)
    return structure

pdbdict = dict()
def GetResidueObj(pdbid, chainid, resname, resnum):
    if pdbid in pdbdict:
        struct = pdbdict[pdbid]
    else:
        struct = GetStructure(pdbid)
        pdbdict[pdbid] = struct
    return GetResidueFromPDB(struct, chainid, resname, resnum)

