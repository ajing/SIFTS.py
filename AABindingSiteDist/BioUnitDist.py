'''
    1. measuring the shortest distance between PDB residues for all biounit files with that PDB
    2. keep all relevant information for that shortest distance
'''

from BioUnitModel import Dist, Residue

BSLineOrder = ["PDBID", "proteinChainID", "resnam", "resnum"]

BSResidueFile = ".txt"
BIOUNITDIR    = "../ligandNet/MOAD2013"

def GetAllPDBDist():
    pdb_bs_residues = AllBSResidueList()
    for eachbiounit in biounitlist:
        bs_residuelist = pdb_bs_residues.getBSResidueList(eachbiounit.split(".")[0])
        for model in biounit_struct:
            for residue in model.get_residues():
                mindist = bs_residuelist.GetCloseDist(residue, biounit_struct)


class AllBSResidueList:
    def __init__(self):
        self.bs_dict = dict()
        self.parserBSFile()

    def parserBSFile(self):
        for line for open(BSResidueFile):
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
