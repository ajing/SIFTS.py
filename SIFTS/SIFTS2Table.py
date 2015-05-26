'''
    SIFTS data to a table, so I can load to MySQL database.
'''

from SIFTSXMLMapControl import processOneXML
from SIFTSXMLMapModel import XMLDIR
import os

## parameters
OUTFILE = "./SIFTS/sift_sql_table_test.txt"
#OUTDIR  = "./SIFTS/sift_sql_test/"
OUTDIR  = "./SIFTS/sift_sql/"
MOADPDB = "./Data/MOADPDB.txt"
MULTIPRO= True

from multiprocessing import Pool
import threading
# thread safe for writing file
mutex_writefile = threading.Lock()


def oneXML2Table(filename):
    #print "filename", filename
    content = []
    protein = processOneXML(filename)
    pdbid   = protein.pdbid
    for eachres in protein.getResidues():
        uniprot = eachres.uniprot
        if uniprot is None:
            continue
        if not all(map(str, [pdbid, eachres.resChain, eachres.resNum, eachres.resName, uniprot.accid, uniprot.resname, uniprot.resnum])):
            continue
        line    = "\t".join(map(str, [pdbid, eachres.resChain, eachres.resNum, eachres.resName, uniprot.accid, uniprot.resname, uniprot.resnum]))
        content.append(line)
    if content:
        fileobj = open(OUTDIR + pdbid, "w")
        fileobj.write("\n".join(content) + "\n")
        fileobj.close()

def FileFilter(inputlist):
    pdblist = []
    #for line in open(OUTFILE):
    #    content = line.strip().split()
    for line in os.listdir(OUTDIR):
        content = line.strip()
        pdb = content
        if len(pdb) != 4:
            print "something wrong with: " + pdb
            continue
        pdblist.append(pdb)
    #print "PDB before set unique: ", len(pdblist)
    pdblist = list(set(pdblist))
    #print "PDB in outdir:", len(pdblist)
    moadlist= []
    for line in open(MOADPDB):
        content  = line.strip()
        moadlist.append(content.lower())
    newlist = []
    moadcopy= list(moadlist)
    #print "moadlist len,",len(moadlist)
    #print "inputlist len,", len(inputlist)
    for eachfile in inputlist:
        pdb = eachfile.split('/')[-1].split('.')[0]
        if not pdb in pdblist:
            newlist.append(eachfile)
        #print "pdb:", pdb
        ## the following code is for BindingMOAD only
        #if pdb in pdblist and pdb in moadlist:
        #    moadcopy.remove(pdb)
        #if not pdb in pdblist and pdb in moadlist:
        #    moadcopy.remove(pdb)
        #    newlist.append(eachfile)
    #print "final list:", len(newlist)
    return newlist

def XML2Table(filedir):
    allcontent = []
    inputlist = []
    #try:
    #    os.remove(OUTFILE)
    #except:
    #    pass
    for root, dirs, files in os.walk(filedir):
        for afile in files:
            afile = os.path.join(root, afile)
            inputlist.append(afile)
    #print "original:", len(inputlist)
    inputlist = FileFilter(inputlist)
    #print len(inputlist)
    if MULTIPRO:
        pool = Pool(processes = 5)
        result = pool.map_async(oneXML2Table, inputlist)
        resulttxt = result.wait()
    else:
        # for single process
        for filename in inputlist:
            oneXML2Table(filename)

if __name__ == "__main__":
    XML2Table(XMLDIR)
