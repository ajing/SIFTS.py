
from SIFTS.SIFTSXMLMapModel import Residue
from GetDistMap import GetSeqSpa
from PDBtools import GetResidueObj

#PDBDIR = "./Data/receptorConCat"
BSDIR  = "./Data/bindingsite2.txt"
OUTDIR = "distaa_biolip_2.txt"
SNPINFO= "./Data/snp_logic_snp.txt"

_getseqspa = GetSeqSpa()

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

    def getminDist(self, tar_resnam, tar_reschain, tar_resnum):
        mindist = float("inf")
        minseqdist = float("inf")
        for each in self.residuelist:
            resname = each.getPDBresName()
            resnum  = each.getPDBresNum()
            try:
                res_obj = GetResidueObj(self._pdbid, self._chainid, resname, resnum)
            except Exception as e:
                print e
                continue
            try:
                spa_dist, seq_dist = _getseqspa.GetSeqSpafor2Residue(self._pdbid, self._chainid, tar_reschain, resname, tar_resnam, resnum, tar_resnum)
            except Exception as e:
                print e
                continue
            if spa_dist < mindist:
                mindist = spa_dist
                minseqdist = seq_dist
        return mindist, minseqdist

    def __str__(self):
        return "PDBID: %s, chainid: %s, bscode: %s, numres: %s;" %(self.PDBID, self.chainid, self.bscode, len(self.residuelist))

