'''
    Get the relation of sequence distance and spatial distance
'''

from SIFTS.SIFTSXMLMapControl import processAllXML, processOneXML, getCorrespondingUniProt
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

class GetSeqSpa:
    def __init__(self):
        self.protein2uniprot = [None, None]
        self.xmldir = XMLDIR

    def GetSeqSpafor2Residue(self, pdbid, chainidA, chainidB, resnamA, resnamB, resnumA, resnumB):
        pdbid   = pdbid.lower()

        try:
            resAstru= GetResidueObj(pdbid, chainidA, resnamA, resnumA)
            resBstru= GetResidueObj(pdbid, chainidB, resnamB, resnumB)
        except Exception as e:
            print e
            return None, None

        if pdbid == self.protein2uniprot[0]:
            protein = self.protein2uniprot[1]
        else:
            protein = self.GetXMLObj(pdbid)
            self.protein2uniprot = [pdbid, protein]

        resAseq = protein.getResidue(chainidA, resnumA)
        resBseq = protein.getResidue(chainidB, resnumB)

        dist_spa = GetSpatialDistance(resAstru, resBstru)
        try:
            dist_seq = resAseq.getSeqDistance(resBseq)
        except:
            dist_seq = False

        #if not dist_seq:
        #    return None, None

        return dist_spa, dist_seq

    def GetXMLObj(self, pdbid):
        filedir = os.path.join(self.xmldir, pdbid[1:3].lower(), pdbid.lower() + ".xml.gz")
        return processOneXML(filedir)

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
                print e
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
    #filedir  = XMLDIR
    #proteins, uniprotdict = processAllXML(filedir)
    #GetSeqSpaAll(uniprotdict)
    seqspa = GetSeqSpa()
    seqspa.GetSeqSpafor2Residue("10gs", "A", "A", "Pro", "Pro", 2, 100)
    seqspa.GetSeqSpafor2Residue("10gs", "A", "A", "", "", 2, 100)
