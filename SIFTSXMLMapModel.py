'''
    Map PDBe residue to SIFTS residue
'''

#PDBDIR = "pdbtest"
#XMLDIR = "xmltest"
PDBDIR = "pdb"
XMLDIR = "xml"

PAIRFILE = "pair.txt"
SAMPLESIZE = 10

class UniProtInfo:
    def __init__(self, accid, resname, resnum):
        self._accid = accid
        self._resname = resname
        try:
            self._resnum  = int(resnum)
        except:
            raise Exception("Cannot parse resnum")

    @property
    def accid(self):
        return self._accid

    @property
    def resname(self):
        return self._resname

    @property
    def resnum(self):
        return self._resnum

    def __eq__(self, uniprot):
        if isinstance(uniprot, UniProtInfo):
            if uniprot.accid == self.accid:
                return True
            else:
                return False
        else:
            return False

    def __ne__(self, uniprot):
        if isinstance(uniprot, UniProtInfo):
            if uniprot.accid != self.accid:
                return True
            else:
                return False
        else:
            return True


class Residue:
    def __init__(self, resnum, resnam, reschain):
        # residue name and number from PDBe
        self.resName = resnam
        self.resNum  = int(resnum)
        self.resChain= reschain
        self.uniprot = None

    def getPDBresName(self):
        return self.resName

    def getPDBresNum(self):
        return self.resNum

    def getPDBresChain(self):
        return self.resChain

    def setUniProtInfo(self, uniprotinfo):
        if not isinstance(uniprotinfo, UniProtInfo):
            raise TypeError("wrong type for uniprot")
        self.uniprot = uniprotinfo

    def getUniProtInfo(self):
        if not self.uniprot is None:
            self.uniprot
        else:
            raise Exception("cannot find uniprot")

    def getSeqDistance(self, res2):
        if self.uniprot == res2.uniprot and not self.uniprot is None:
            return abs(self.uniprot.resnum - res2.uniprot.resnum)
        else:
            return False


    def __str__(self):
        return "%s,%s;" %(self.resName, self.resNum)

class Protein:
    def __init__(self, pdbid):
        self.resList = []
        self.uniprots= []
        self.pdbid = pdbid

    def getProteinID(self):
        return self.pdbid

    def appendNewResidue(self, residue):
        if not isinstance(residue, Residue):
            raise TypeError("wrong type for residue")
        self.resList.append(residue)

    def appendNewUniProt(self, uniprot):
        if not isinstance(uniprot, UniProtInfo):
            raise TypeError("wrong type for uniprot")
        if not uniprot in self.uniprots:
            self.uniprots.append(uniprot)

    def getResidues(self):
        return self.resList

    def getUniProts(self):
        return self.uniprots

    def length(self):
        return len(self.resList)

    def iterateResidue(self):
        return iter(self.resList)

    def __str__(self):
        prints = self.pdbid + "\t"
        for each in self.resList:
            prints += str(each)
        return prints
