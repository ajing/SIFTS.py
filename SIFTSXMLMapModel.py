'''
    Map PDBe residue to SIFTS residue
'''

class UniProtInfo:
    def __init__(self, accid, resname, resnum):
        self.accid = accid
        self.resname = resname
        self.resnum  = resnum

    def getAccid():
        return self.accid

    def getResName():
        return self.resname

    def getResNum():
        return self.resnum

    def __cmp__(self, uniprot):
        if isinstance(uniprot, UniProtInfo):
            if uniprot.getAccid() == self.accid:
                return True
            else:
                return False
        else:
            return False


class Residue:
    def __init__(self, resnum, resnam, reschain):
        # residue name and number from PDBe
        self.resName = resnam
        self.resNum  = resnum
        self.resChain= reschain

    def getPDBresName(self):
        return self.resName

    def getPDBresNum(self):
        return self.resNum

    def getPDBresChain(self):
        return self.resChain

    def setUniProtInfo(self, accid, resnum, resname):
        self.uniprot = UniProtInfo(accid, resname, resnum)

    def getUniProtInfo(self):
        try:
            self.uniprot
        except:
            raise Exception("cannot find uniprot")

    def getSeqDistance(self, res2):
        if self.uniprot == res2.uniprot:
            return self.uniprot.getResNum() - res2.uniprot.getResNum()
        else:
            return False


    def __str__(self):
        return "%s,%s;" %(self.resName, self.resNum)

class Protein:
    def __init__(self, pdbid):
        self.resList = []
        self.pdbid = pdbid

    def getProteinID(self):
        return self.pdbid

    def appendNewResidue(self, residue):
        if not isinstance(residue, Residue):
            raise TypeError("wrong type for residue")
        self.resList.append(residue)

    def getResidues():
        return self.resList

    def iterateResidue(self):
        return iter(self.resList)

    def __str__(self):
        prints = self.pdbid + "\t"
        for each in self.resList:
            prints += str(each)
        return prints
