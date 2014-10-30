'''
    Get the relation of sequence distance and spatial distance
'''

from SIFTSXMLMapControl import processAllXML
from PDBtools import GetResidueObj
from SIFTSXMLMapModel import XMLDIR

def GetSpatialDistance(res1, res2):
    return res1['CA'] - res2['CA']

def GetSeqSpaPair(protein):
    reslist  = protein.getResidues()
    pairlist = []
    pdb      = protein.getProteinID()
    reslen   = len(reslist)
    for i in range(reslen - 1):
        res1    = reslist[i]
        chainid = res1.getPDBresChain()
        resnum  = res1.getPDBresNum()
        resnam  = res1.getPDBresName()
        try:
            res1stru= GetResidueObj(pdb, chainid, resnam, resnum)
        except Exception as e:
            #print e
            continue
        for j in range(i + 1, reslen):
            res2    = reslist[j]
            chainid = res2.getPDBresChain()
            resnum  = res2.getPDBresNum()
            resnam  = res2.getPDBresName()
            try:
                res2stru= GetResidueObj(pdb, chainid, resnam, resnum)
            except Exception as e:
                #print e
                continue
            dist = res1.getSeqDistance(res2)
            if dist:
                pairlist.append((GetSpatialDistance(res1stru, res2stru), dist))
    return pairlist

def GetSeqSpaAll(uniprotdict):
    pairlistall = []
    for protein in uniprotdict.values():
        pairlistall += GetSeqSpaPair(protein)
    return pairlistall


def SavePairList(pairlist):
    fileobj = open("pair.txt", "w")
    for each in pairlist:
        print each
        spatialdist, residuedist = each
        line = "%s\t%s\n" % (spatialdist, residuedist)
        fileobj.write(line)
    fileobj.close()

if __name__ == "__main__":
    filedir  = "xml"
    proteins, uniprotdict = processAllXML(filedir)
    pairlist = GetSeqSpaAll(uniprotdict)
    SavePairList(pairlist)
