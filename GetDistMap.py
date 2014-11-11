'''
    Get the relation of sequence distance and spatial distance
'''

from SIFTS.SIFTSXMLMapControl import processAllXML
from PDBtools import GetResidueObj
from SIFTS.SIFTSXMLMapModel import XMLDIR, PAIRFILE
import os

def GetSpatialDistance(res1, res2):
    try:
        return res1['CA'] - res2['CA']
    except:
        print res1
        print res2
        raise Exception("cannot find CA")

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
            print e
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
            dist_seq = res1.getSeqDistance(res2)
            try:
                dist_spa = GetSpatialDistance(res1stru, res2stru)
            except:
                continue
            if dist_seq and dist_spa:
                pairlist.append((dist_spa, dist_seq))
    return pairlist

def GetSeqSpaAll(uniprotdict):
    try:
        os.remove(PAIRFILE)
    except:
        pass
    for protein in uniprotdict.values():
        pairlist = GetSeqSpaPair(protein)
        SavePairList(pairlist)

def SavePairList(pairlist):
    fileobj = open(PAIRFILE, "a")
    for each in pairlist:
        spatialdist, residuedist = each
        line = "%s\t%s\n" % (spatialdist, residuedist)
        fileobj.write(line)
    fileobj.close()

if __name__ == "__main__":
    filedir  = XMLDIR
    proteins, uniprotdict = processAllXML(filedir)
    GetSeqSpaAll(uniprotdict)
