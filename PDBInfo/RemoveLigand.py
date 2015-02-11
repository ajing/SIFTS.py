'''
    Remove ligand in biounit file
'''

from Bio.PDB import PDBParser, PDBIO
import os

LIGANDLIST = "../Data/ligandlist.txt"
BIOUNITDIR = "../../ligandNet/2013_biounits"
BIOSTRDIR  = "../../ligandNet/2013_biounits_noligand"
LIGANDCOL  = ["PDBID", "LigName", "ChainID", "ResNum", "Validity", "LigInfo"]

def ProcessLigFile(infile):
    pdbligand = dict()
    for line in open(infile):
        content = line.strip().split("\t")
        contdict = dict(zip(LIGANDCOL, content))
        pdb      = contdict["PDBID"]
        if pdb in pdbligand:
            pdbligand[pdb].append(contdict)
        else:
            pdbligand[pdb] = [contdict]
    return pdbligand


def LongLigand(chain, residue, longligand):
    lig_len = len(longligand["LigName"].split())
    if chain.id == longligand["ChainID"] and residue.id[1] - int(longligand["ResNum"]) < lig_len * 1.1 and residue.resname in longligand["LigName"]:
        chain.detach_child(residue.id)

def RemoveLigandsOneBioUnit(biounit, ligandlist):
    # ligandlist is a residue list with residue chain id, name and residue number
    p = PDBParser(PERMISSIVE = 1)
    pdbname= biounit.split("/")[-1]
    try:
        models = p.get_structure(pdbname, biounit)
    except:
        return None
    print ligandlist
    #for model in models:
    #    for chain in model:
    #        for residue in chain:
    #            print residue
    for rligand in ligandlist:
        for model in models:
            for chain in model:
                for residue in list(chain):
                    if chain.id == rligand["ChainID"] and int(rligand["ResNum"]) == residue.id[1]:
                        chain.detach_child(residue.id)
                    elif residue.id[0] == "W":
                        chain.detach_child(residue.id)
                    elif len(rligand["LigName"].split()) > 1 and int(rligand["ResNum"]) <= residue.id[1]:
                        LongLigand(chain, residue, rligand)
    io = PDBIO()
    io.set_structure(models)
    filepath = os.path.join(BIOSTRDIR, models.id)
    io.save(filepath)

def RemoveLigandsAll():
    pdbligand = ProcessLigFile(LIGANDLIST)
    for eachbiounit in os.listdir(BIOUNITDIR):
        biounit = os.path.join(BIOUNITDIR, eachbiounit)
        pdb     = eachbiounit.split(".")[0].upper()
        RemoveLigandsOneBioUnit(biounit, pdbligand[pdb])

if __name__ == "__main__":
    #ProcessLigFile(LIGANDLIST)
    RemoveLigandsAll()
