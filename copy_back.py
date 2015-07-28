'''
    Copy the file back to pdb folder
'''

import shutil
import os

ORIG_DIR = "pdb_gz"
DEST_DIR = "pdb"


def CopyFile(filename):
    if not filename.endswith(".ent.gz"):
        return
    srcfile  = os.path.join(ORIG_DIR, filename)
    dstdir   = os.path.join(DEST_DIR, filename[4:6], filename)
    #print srcfile
    #print dstdir
    shutil.copy(srcfile, dstdir)

def main():
    filelist = os.listdir(ORIG_DIR)
    for afile in filelist:
        CopyFile(afile)

if __name__ == "__main__":
    CopyFile("pdb10gs.ent.gz")
    main()
