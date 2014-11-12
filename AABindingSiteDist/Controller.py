'''
    Calculate the distance between binding site to other residues
    python -m AABindingSiteDist.Controller
'''

from PDBtools import GetResidueObj
from Model import BindingSite
from SIFTS.SIFTSXMLMapModel import Residue
from Bio.PDB.Polypeptide import one_to_three

BSLineOrder = ["PDBID", "ChainID", "BSID", "LigName", "LigChain", "BSRes"]
SNPLineOrder = ["PDBID", "ChainID", "SNPID", "ResNum", "ResName"]

def GetAllDist(bslist):
    file_obj = open("distaa.txt", "w")
    for eachbs in bslist:
        pdbid    = eachbs.pdbid
        chainid  = eachbs.chainid
        bscode   = eachbs.bscode
        reslist  = GetResidueObj(pdbid)
        print pdbid, chainid, bscode
        for residue in reslist:
            resname = residue.get_resname()
            reschain= residue.get_full_id()[2]
            resnum  = residue.id[1]
            mindist = eachbs.getminDist(residue)
            info    = [pdbid, chainid, bscode, resname, reschain, resnum, mindist]
            line    = "\t".join(map(str, info))
            file_obj.write(line + "\n")
    file_obj.close()

def BSParser(infile):
    bslist = []
    for line in open(infile):
        content = line.strip().split('\t')
        pdbid   = content[BSLineOrder.index("PDBID")].lower()
        chainid = content[BSLineOrder.index("ChainID")]
        bscode  = content[BSLineOrder.index("BSID")]
        ligname = content[BSLineOrder.index("LigName")]
        ligchain= content[BSLineOrder.index("LigName")]
        bsres   = content[BSLineOrder.index("BSRes")]
        newbs   = BindingSite(pdbid, chainid, bscode, ligchain, ligname)
        for eachres in bsres.split():
            resname = one_to_three(eachres[0])
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
    return bslist

def SNPResParser(infile):
    snpres = dict()
    for line in open(infile):
        pdbid   = content[BSLineOrder.index("PDBID")].lower()
        chainid = content[BSLineOrder.index("ChainID")]
        snpid   = content[BSLineOrder.index("SNPID")]
        resnum  = content[BSLineOrder.index("ResNum")]
        resnam  = content[BSLineOrder.index("ResName")]
        snpres[(pdbid, chainid)] = (snpid, resnum, resnam)
    return snpres

'''
class FilterSNP:
    def __init__(self, snpres):
        self.snpdict = snpres

    def isSNP(self, pdbid, chainid, resnum):
        if not (pdbid, chainid) in
        snpid, resnum, resnam = snpres[(pdbid)]
'''


if __name__ == "__main__":
    bslist = BSParser("./Data/bindingsite.txt")
    GetAllDist(bslist)
