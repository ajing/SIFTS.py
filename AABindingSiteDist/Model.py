
from SIFTS.SIFTSXMLMapModel import Residue
from GetDistMap import GetSpatialDistance
from PDBtools import GetResidueObj

PDBDIR = "./Data/receptorConCat"
BSDIR  = "./Data/bindingsite2.txt"
OUTDIR = "distaa_biolip.txt"

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
        if not isinstance(res, Residue):
            raise Exception("this is not a residue object in SIFT model")
        self.residuelist.append(res)
        resname = res.getPDBresName()
        resnum = res.getPDBresNum()

    def getminDist(self, res):
        mindist = float("inf")
        for each in self.residuelist:
            resname = each.getPDBresName()
            resnum  = each.getPDBresNum()
            try:
                res_obj = GetResidueObj(self._pdbid, self._chainid, resname, resnum)
            except Exception as e:
                print e
                continue
            try:
                dist = GetSpatialDistance(res_obj, res)
            except Exception as e:
                print e
                continue
            if dist < mindist:
                mindist = dist
        return mindist

    def __str__(self):
        return "PDBID: %s, chainid: %s, bscode: %s, numres: %s;" %(self.PDBID, self.chainid, self.bscode, len(self.residuelist))

