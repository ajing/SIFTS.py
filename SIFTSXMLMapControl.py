'''
    Save data to MapObject
'''

from xml.dom import minidom
import ntpath
from SIFTSXMLMapModel import Protein, Residue

def processOneXML(filename):
    xmldoc    = minidom.parse(filename)
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
                    uniprotid = cref.getAttribute('dbAccessionId')
                    uniprotnum = cref.getAttribute('dbResNum')
                    uniprotnam = cref.getAttribute('dbResName')
                    residue.setUniProtInfo(uniprotid, uniprotnum, uniprotnam)
                    break
    return protein

def processAllXML(filedir):
    proteins = dict()
    for file in os.listdir(filedir):
        protein = processOneXML(file)
        proteins[protein.getProteinID] = protein
    return proteins

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
