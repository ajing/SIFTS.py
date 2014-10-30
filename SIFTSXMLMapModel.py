'''
    Map PDBe residue to SIFTS residue
'''

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
        self.uniprot = {"accid": accid, "residueName:": resname, "ResidueNum:": resnum}

    def getUniProtInfo(self):
        try:
            self.uniprot
        except:
            raise Exception("cannot find uniprot")

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

    def length(self):
        return len(self.resList)

    def iterateResidue(self):
        return iter(self.resList)

    def __str__(self):
        prints = self.pdbid + "\t"
        for each in self.resList:
            prints += str(each)
        return prints
