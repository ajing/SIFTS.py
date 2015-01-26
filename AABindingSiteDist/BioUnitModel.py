'''
    Model for BioUnitDist.py
'''

class Residue:
    def __init__(self, pdbid, chainid, resnum, resnam):
        self._pdbid   = pdbid
        self._chainid = chainid
        self._resnum  = resnum
        self._resnam  = resnam

    @property
    def pdbid(self):
        return self._pdbid

    @property
    def chainid(self):
        return self._chainid

    @property
    def resnum(self):
        return self._resnum

    @property
    def resnam(self):
        return self._resnam

class Dist:
    def __init__(self, biounit, chain1, chain2, resnum1, resnum2, resnam1, resnam2):
        self._biounit = None
        self._chainid1 = None
        self._chainid2 = None
        self._resnum1  = None
        self._resnum2  = None
        self._resnam1  = None
        self._resnam2  = None

    @property
    def biounit(self):
        return self._biounit

    @property
    def chainid1(self):
        return self._chainid1

    @property
    def chainid2(self):
        return self._chainid2

    @property
    def resnum1(self):
        return self._resnum1

    @property
    def resnum2(self):
        return self._resnum2

    @property
    def resnam1(self):
        return self._resnam1

    @property
    def resnam2(self):
        return self._resnam2

def GetAllBioUnitDist(pdbid, res1, res2):
    for eachbiounit in GetAllBioUnit(pdbid):
        dist = GetResidueDist(eachbiounit, res1, res2)
        line = [eachbiounit, chainid1, resnam1, resnum1, chainid2, resnam2, resnum2, dist]

def GetResidueFromBioUnit(biounit_str, chainid, resname, resnum):
    for model in biounit_str:
        model_resnam = model[chainid][int(resnum)].get_resname()
        if model_resnam != resname:
            print "in %s, %s doesn't match %s." % (biounit_str.get_id(), resname, model_resnam)
        yield model.get_id(), model[chainid][int(resnum)]


class BSResidueList
    def __init__(self, pdbid):
        self.pdb = pdbid
        self.bs_res_list = []

    def addBSResidue(self, residue):
        self.bs_res_list.append(residue)

    def getCloseDist(self, residue, biounit_struct):
        mindist = None
        for bs_res in self.bs_res_list:
            model_id, bs_res_obj = GetResidueFromBioUnit(biounit_struct, bs_res.chainid, bs_res.resnam, bs_res.resnum)
            distance   = bs_res_obj - residue
            if mindist > distance or mindist is None:
                    mindist = distance
        return mindist
