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
    proteins = []
    for file in os.listdir(filedir):
        proteins.append(processOneXML(file))

if __name__ == "__main__":
    processOneXML("200l.xml")
