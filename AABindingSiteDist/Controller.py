'''
    Calculate the distance between binding site to other residues
    python -m AABindingSiteDist.Controller
'''

from PDBtools import GetResidueObj, GetFilewithPDB
from Model import BindingSite, BSDIR, OUTDIR
from SIFTS.SIFTSXMLMapModel import Residue
from Bio.PDB.Polypeptide import one_to_three
import os

BSLineOrder = ["PDBID", "ChainID", "BSID", "LigName", "LigChain", "BSRes"]

def GetAllDist(bslist):
    file_obj = open(OUTDIR, "w")
    for eachbs in bslist:
        pdbid    = eachbs.pdbid
        chainid  = eachbs.chainid
        bscode   = eachbs.bscode
        try:
            reslist  = GetResidueObj(pdbid)
        except Exception as e:
            print e
            print "cannot find reslist for " + pdbid
            continue

        print pdbid, chainid, bscode
        for residue in reslist:
            resname = residue.get_resname()
            reschain= residue.get_full_id()[2]
            resnum  = residue.id[1]
            min_spa_dist, min_seq_dist = eachbs.getminDist(resname, reschain, resnum)
            info    = [pdbid, chainid, bscode, resname, reschain, resnum, min_spa_dist, min_seq_dist]
            line    = "\t".join(map(str, info))
            file_obj.write(line + "\n")
    file_obj.close()

def PDBFileExist(pdbid):
    filedir = GetFilewithPDB(pdbid)
    if os.path.isfile(filedir):
        return True
    else:
        return False

def BSParser(infile):
    bslist = []
    for line in open(infile):
        content = line.strip().split('\t')
        pdbid   = content[BSLineOrder.index("PDBID")].lower()
        chainid = content[BSLineOrder.index("ChainID")]
        bscode  = content[BSLineOrder.index("BSID")]
        ligname = content[BSLineOrder.index("LigName")]
        ligchain= content[BSLineOrder.index("LigChain")]
        bsres   = content[BSLineOrder.index("BSRes")]
        newbs   = BindingSite(pdbid, chainid, bscode, ligchain, ligname)
        for eachres in bsres.split():
            try:
                resname = one_to_three(eachres[0])
            except:
                print "wrong bindingsite res: " + eachres
                continue
            try:
                resnum = int(eachres[1:])
            except:
                continue
                #raise Exception("convert %s to number" % eachres[1:])
            residue = Residue(resnum, resname, chainid)
            try:
                newbs.appendResidue(residue)
            except Exception as e:
                print e
                continue
        bslist.append(newbs)

        if not PDBFileExist(pdbid):
            print "Cannot find file " + pdbid
    return bslist

def PrintBSList(bslist):
    for eachbs in bslist:
        for eachres in eachbs.residuelist:
            content = [eachbs.pdbid, eachbs._ligchain, eachbs._ligname, eachbs.bscode, eachres.getPDBresChain(), eachres.getPDBresName(), eachres.getPDBresNum()]
            print "\t".join(map(str, content))

if __name__ == "__main__":
    bslist = BSParser(BSDIR)
    PrintBSList(bslist)
    #GetAllDist(bslist)
