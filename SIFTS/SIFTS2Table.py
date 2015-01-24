'''
    SIFTS data to a table, so I can load to MySQL database.
'''

from SIFTSXMLMapControl import processOneXML
from SIFTSXMLMapModel import XMLDIR
import os

## parameters
OUTFILE = "./SIFTS/sift_sql_table.txt"

from multiprocessing import Pool
import threading
# thread safe for writing file
mutex_writefile = threading.Lock()

def oneXML2Table(filename):
    print "filename", filename
    content = []
    protein = processOneXML(filename)
    pdbid   = protein.pdbid
    for eachres in protein.getResidues():
        uniprot = eachres.uniprot
        if uniprot is None:
            continue
        line    = "\t".join(map(str, [pdbid, eachres.resChain, eachres.resNum, eachres.resName, uniprot.accid, uniprot.resname, uniprot.resnum]))
        content.append(line)
    if content:
        mutex_writefile.acquire()
        fileobj = open(OUTFILE, "a")
        fileobj.write("\n".join(content) + "\n")
        fileobj.close()
        mutex_writefile.release()

def XML2Table(filedir):
    allcontent = []
    inputlist = []
    try:
        os.remove(OUTFILE)
    except:
        pass
    for root, dirs, files in os.walk(filedir):
        for afile in files:
            afile = os.path.join(root, afile)
            inputlist.append(afile)
    pool = Pool(processes = 7)
    result = pool.map_async(oneXML2Table, inputlist)
    resulttxt = result.wait()

if __name__ == "__main__":
    XML2Table(XMLDIR)
