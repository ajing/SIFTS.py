'''
    SIFTS data to a table, so I can load to MySQL database.
'''

from SIFTSXMLMapControl import processOneXML
from SIFTSXMLMapModel import XMLDIR
import os

## parameters
OUTFILE = "./SIFTS/sift_sql_table.txt"

def oneXML2Table(filename, fileobj):
    content = []
    protein = processOneXML(filename)
    pdbid   = protein.pdbid
    print pdbid
    for eachres in protein.getResidues():
        uniprot = eachres.uniprot
        if uniprot is None:
            continue
        line    = "\t".join(map(str, [pdbid, eachres.resChain, eachres.resNum, eachres.resName, uniprot.accid, uniprot.resname, uniprot.resnum]))
        content.append(line)
    fileobj.write("\n".join(content) + "\n")

def XML2Table(filedir):
    allcontent = []
    outobj = open(OUTFILE, "w")
    for root, dirs, files in os.walk(filedir):
        for afile in files:
            afile = os.path.join(root, afile)
            oneXML2Table(afile, outobj)
    outobj.close()

if __name__ == "__main__":
    XML2Table(XMLDIR)
