'''
    Save data to MapObject
'''

from xml.dom import minidom
import ntpath
import os
from SIFTSXMLMapModel import UniProtInfo, Protein, Residue
from common import length
import gzip

def processOneXML(filename):
    xmldoc    = minidom.parse(gzip.open(filename))
    entitylist= xmldoc.getElementsByTagName('entity')
    pdb       = ntpath.basename(filename)[:4]
    protein   = Protein(pdb)
    for entity in entitylist:
        residuelist = entity.getElementsByTagName('residue')
        chainid     = entity.getAttribute('entityId')
        for residue in residuelist:
            resnum    = residue.getAttribute('dbResNum')
            resname   = residue.getAttribute('dbResName')
            crossrefs = residue.getElementsByTagName('crossRefDb')
            residue   = Residue(resnum, resname, chainid)
            protein.appendNewResidue(residue)
            for cref in crossrefs:
                if cref.getAttribute('dbSource') == "UniProt":
                    uniprotid  = cref.getAttribute('dbAccessionId')
                    uniprotnum = cref.getAttribute('dbResNum')
                    uniprotnam = cref.getAttribute('dbResName')
                    uniprotinfo = UniProtInfo(uniprotid, uniprotnam, uniprotnum)
                    residue.setUniProtInfo(uniprotinfo)
                    protein.appendNewUniProt(uniprotinfo)
                    break
    return protein

def processAllXML(filedir):
    proteins = dict()
    uniprot  = dict()
    for root, dirs, files in os.walk(filedir):
        for afile in files:
            afile = os.path.join(root, afile)
            protein = processOneXML(afile)
            proteins[protein.getProteinID()] = protein
            uniprotinfos = protein.getUniProts()
            for uniprotinfo in uniprotinfos:
                accid       = uniprotinfo.accid
                if accid in uniprot:
                    if protein.length() > uniprot[accid].length():
                        uniprot[accid] = protein
                else:
                    uniprot[accid] = protein
    return (proteins, uniprot)

def getCorrespondingUniProt(proteins, pdbid, chainid, resnum):
    try:
        protein = proteins[pdbid]
    except:
        raise Exception("Cannot find corresponding UniProt for:" + pdbid)
    for res in protein:
        if chainid == res.getPDBresChain() and resnum == res.getPDBresNum():
            return res.getUniProtInfo(pdbid, chainid, resnum)

    raise Exception("Cannot find corresponding UniProt for: %s, %s, %s" % (pdbid, chainid, resnum))


if __name__ == "__main__":
    processOneXML("200l.xml")
