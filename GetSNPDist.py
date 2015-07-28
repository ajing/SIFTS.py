'''
    Get the spatial and sequence distance between SNP residues.
'''

from GetDistMap import GetSeqSpa
from SIFTS.SIFTSXMLMapModel import Residue

COLNAME = ["pdbid", "snpid", "chainid", "resnum", "resnam", "alt_residue", "chr", "strand", "chr_position", "secondary", "solubility", "ligandcode", "conservation", "domain_interface", "genename", "SwissProt_AC", "FTID", "aachange", "vartype", "dbSNPid", "disease_name"]

SNP_INFILE = "./Data/snpannotate.txt"
OUTFILE    = "./Result/snp_dist_test.txt"

def SNP_pair_dist(pdbid, reslist):
    getsesp = GetSeqSpa()
    # residues is a list of residue information
    reslen  = len(reslist)
    content = ""
    for i in range(reslen - 1):
        res1    = reslist[i]
        chainid1 = res1.getPDBresChain()
        resnum1  = res1.getPDBresNum()
        resnam1  = res1.getPDBresName()
        snpid1   = res1.snpid
        for j in range(i + 1, reslen):
            res2    = reslist[j]
            chainid2 = res2.getPDBresChain()
            resnum2  = res2.getPDBresNum()
            resnam2  = res2.getPDBresName()
            snpid2   = res2.snpid
            dist_spa, dist_seq = getsesp.GetSeqSpafor2Residue(pdbid, chainid1, chainid2, resnam1, resnam2, resnum1, resnum2)
            if dist_spa is None:
                continue
            line = [pdbid, snpid1, chainid1, resnam1, resnum1, snpid2, chainid2, resnam2, resnum2, dist_spa, dist_seq]
            content += "\t".join(map(str, line)) + "\n"
    return content

def Infile2DictList(infile):
    dictlist = []
    #count  = 0
    for line in open(infile):
        content = line.strip().split("\t")
        dictlist.append(dict(zip(COLNAME, content)))
        #count += 1
        #if count > 500:
        #    break
    return dictlist

def ParseInput(snploc_file):
    dictlist = Infile2DictList(snploc_file)
    pdbdict  = dict()
    for oneline in dictlist:
        pdbid   = oneline["pdbid"]
        residue = Residue(int(oneline["resnum"]), oneline["resnam"], oneline["chainid"], oneline["snpid"])
        if pdbid in pdbdict:
            pdbdict[pdbid].append(residue)
        else:
            pdbdict[pdbid] = [residue]
    return pdbdict


def SNP_pair_all(pdbdict):
    outobj = open(OUTFILE, "w")
    #for pdbid in pdbdict:
    #    content = SNP_pair_dist(pdbid, pdbdict[pdbid])
    #    outobj.write(content)
    print pdbdict.keys()
    content = SNP_pair_dist("1HAQ", pdbdict["1HAQ"])
    outobj.close()

def main():
    pdbdict = ParseInput(SNP_INFILE)
    SNP_pair_all(pdbdict)

if __name__ == "__main__":
    main()
