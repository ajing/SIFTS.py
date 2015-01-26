'''
    1. measuring the shortest distance between PDB residues for all biounit files with that PDB
    2. keep all relevant information for that shortest distance
'''

from BioUnitModel import Dist, Residue, BSResidueList
from Bio.PDB import PDBParser
import os

from multiprocessing import Pool
import threading
# thread safe for write file
mutex_file = threading.Lock()

BSLineOrder = ["PDBID", "proteinChainID", "resnam", "resnum"]

BSResidueFile = "./Data/bs_list.txt"
BIOUNITDIR    = "../ligandNet"
OUTFILE       = "./Result/bs_dist_biounit.txt"
#BIOUNITDIR    = "../ligandNet/2013_biounits"

def GetAllPDBDist():
    pdb_bs_residues = AllBSResidueList()
    biounitlist     = os.listdir(BIOUNITDIR)
    GetOneBioUnitDist.bs_residues = pdb_bs_residues
    try:
        os.remove(OUTFILE)
    except:
        pass
    pool = Pool(processes = 5)
    result = pool.map_async(GetOneBioUnitDist, biounitlist)
    result.wait()

def GetOneBioUnitDist(biounit):
    bs_residuelist = GetOneBioUnitDist.bs_residues.getBSResidueList(biounit.split(".")[0])
    biounit_struct = GetStructure(biounit)
    content = []
    for model in biounit_struct:
        for residue in model.get_residues():
            mindist = bs_residuelist.getCloseDist(residue, biounit_struct)
            if not mindist is None:
                mindist.model1 = model.get_id()
                content.append(mindist.printInfo())
    mutex_file.acquire()
    outobj = open(OUTFILE, "a")
    outobj.write("\n".join(content))
    outobj.close()
    mutex_file.release()

def GetStructure(biounit):
    parser = PDBParser(PERMISSIVE = 1)
    pdbobj = open(os.path.join(BIOUNITDIR, biounit))
    structure = parser.get_structure(biounit, pdbobj)
    return structure

class AllBSResidueList:
    def __init__(self):
        self.bs_dict = dict()
        self.parserBSFile()

    def parserBSFile(self):
        for line in open(BSResidueFile):
            content = line.strip().split("\t")
            pdbid   = content[BSLineOrder.index("PDBID")]
            prochain= content[BSLineOrder.index("proteinChainID")]
            resnam  = content[BSLineOrder.index("resnam")]
            resnum  = content[BSLineOrder.index("resnum")]
            residue = Residue(pdbid, prochain, resnum, resnam)
            if pdbid in self.bs_dict:
                self.bs_dict[pdbid].addBSResidue(residue)
            else:
                self.bs_dict[pdbid] = BSResidueList(pdbid)
                self.bs_dict[pdbid].addBSResidue(residue)

    def getBSResidueList(self, pdbid):
        return self.bs_dict[pdbid]

def main():
    GetAllPDBDist()

if __name__ == "__main__":
    main()
