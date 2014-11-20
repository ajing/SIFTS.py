'''
    Concatenate files into one PDB
'''
import os

OUTDIR = "./receptorConCat"
INDIR = "./receptor"

def GetUniquePDB(flist):
    pdblist = []
    for fname in flist:
        pdb = fname[:4]
        pdblist.append(pdb)
    return list(set(pdblist))

def GetFileList(pdbid, filelist):
    flist = []
    for fname in filelist:
        if fname.startswith(pdbid):
            flist.append(os.path.join(INDIR, fname))
    return flist

def ConcateFiles(filelist):
    finalfile = os.path.split(filelist[0])[-1][:4] + ".pdb"
    print finalfile
    filepath  = os.path.join(OUTDIR, finalfile)
    with open(filepath, 'w') as outfile:
        for fname in filelist:
            with open(fname) as infile:
                for line in infile:
                    if not line.startswith("TER"):
                        outfile.write(line)
        outfile.write("TER\n")

def Concate():
    filelist  = os.listdir(INDIR)
    uniquepdb = GetUniquePDB(filelist)
    for eachpdb in uniquepdb:
        flist = GetFileList(eachpdb, filelist)
        ConcateFiles(flist)

if __name__ == "__main__":
    Concate()
