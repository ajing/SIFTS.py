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

    def __str__(self):
        return "%s, %s, %s, %s" %(self._pdbid,self._chainid, self._resnum, self._resnam)

class Dist(object):
    def __init__(self, biounit, chain1, chain2, resnum1, resnum2, resnam1, resnam2, distance):
        self._biounit = biounit
        self._model1  = None
        self._model2  = None
        self._chainid1 = chain1
        self._chainid2 = chain2
        self._resnum1  = resnum1
        self._resnum2  = resnum2
        self._resnam1  = resnam1
        self._resnam2  = resnam2
        self._distance = distance

    @property
    def biounit(self):
        return self._biounit

    @property
    def model1(self):
        return self._model1

    @model1.setter
    def model1(self, modelnum):
        self._model1 = modelnum

    @property
    def model2(self):
        return self._model2

    @model2.setter
    def model2(self, modelnum):
        self._model2 = modelnum

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

    @property
    def distance(self):
        return self._distance

    def __str__(self):
        return ",".join("%s: %s   " % item for item in self.__dict__.items())

    def printInfo(self):
        return "\t".join("%s" % item[1] for item in self.__dict__.items())


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

def GetDistance(res1, res2):
    try:
        return res1['CA'] - res2['CA']
    except:
        return None

class BSResidueList:
    def __init__(self, pdbid):
        self.pdb = pdbid
        self.bs_res_list = []

    def addBSResidue(self, residue):
        self.bs_res_list.append(residue)

    def getCloseDist(self, residue, biounit_struct):
        mindist = None
        for bs_res in self.bs_res_list:
            for model_id, bs_res_obj in GetResidueFromBioUnit(biounit_struct, bs_res.chainid, bs_res.resnam, bs_res.resnum):
                distance   = GetDistance(bs_res_obj, residue)
                dist       = Dist(biounit_struct.get_id(), residue.get_parent().get_id(), bs_res.chainid, residue.get_id()[1], bs_res.resnum, residue.get_resname(), bs_res.resnam, distance)
                dist.model2 = model_id
                if distance is None:
                    continue
                if mindist is None or mindist.distance > distance:
                    mindist = dist
        return mindist
