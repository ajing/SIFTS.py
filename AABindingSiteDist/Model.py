
from SIFTS.SIFTSXMLMapModel import Residue
from GetDistMap import GetSpatialDistance

class BindingSite:
    def __init__(self, pdbid, chainid, bscode, ligchainid, ligname):
        self._pdbid   = pdbid
        self._chainid = chainid
        self._bscode  = bscode
        self._ligchain= ligchainid
        self._ligname = ligname
        #self._ligse   = ligserial
        self.residuelist = []

    @property
    def pdbid(self):
        return self._pdbid

    @property
    def chainid(self):
        return self._chainid

    @property
    def bscode(self):
        return self._bscode

    def appendResidue(self, res):
        self.residuelist.append(res)

    def getminDist(self, res):
        mindist = float("inf")
        for each in self.residuelist:
            dist = GetSpatialDistance(each, res)
            if dist < mindist:
                mindist = dist
        return mindist

    def __str__(self):
        return "PDBID: %s, chainid: %s, bscode: %s, numres: %s;" %(self.PDBID, self.chainid, self.bscode, len(self.residuelist))

