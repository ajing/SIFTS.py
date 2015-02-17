'''
    Filter all distance, only keep mutated residue
'''

from Model import OUTDIR, SNPINFO

SNPLineOrder = ["PDBID", "ChainID", "SNPID", "ResNum", "ResName"]
DISTOrder = ["PDBID", "ChainID", "BSID", "ResName", "ResChain", "ResNum", "Distance", "ResDist"]

#from SCOPData import protein_letters_3to1

def SNPResParser(infile):
    snpres = dict()
    for line in open(infile):
        content = line.strip().split()
        pdbid   = content[SNPLineOrder.index("PDBID")].lower()
        chainid = content[SNPLineOrder.index("ChainID")]
        snpid   = content[SNPLineOrder.index("SNPID")]
        resnum  = content[SNPLineOrder.index("ResNum")]
        resnam  = content[SNPLineOrder.index("ResName")]
        snpres[(pdbid, chainid, resnum)] = (snpid, resnam)
    return snpres

class FilterSNP:
    def __init__(self, snpres):
        self.snpdict = snpres

    def isSNP(self, pdbid, chainid, resnum):
        if not (pdbid, chainid, resnum) in self.snpdict:
            return False
        snpid, snp_resnam = self.snpdict[(pdbid, chainid, resnum)]
        return snpid

def filterFile(infile, filtersnp):
    outfile = infile + "_filtered"
    outobj  = open(outfile, "w")
    for line in open(infile):
        content = line.strip().split()
        pdbid   = content[DISTOrder.index("PDBID")].lower()
        reschain= content[DISTOrder.index("ResChain")]
        resnum  = content[DISTOrder.index("ResNum")]
        snpid   =  filtersnp.isSNP(pdbid, reschain, resnum)
        distance= content[DISTOrder.index("Distance")]
        if snpid and not distance == "inf":
            newline = line.strip() + "\t" + snpid + "\n"
            outobj.write(newline)
    outobj.close

if __name__ == "__main__":
    snpdict   = SNPResParser(SNPINFO)
    filterobj = FilterSNP(snpdict)
    #filterFile("../distaa.txt", filterobj)
    filterFile(OUTDIR, filterobj)
