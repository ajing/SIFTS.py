'''
  Get the nucleotide sequence for each protein
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import seq3

UNIPROT_ID = "uniprot_protein.txt"

def GetTranslatedSeq(seq):
    seq_s = seq.strip().replace(" ", "").replace("\n", "")
    ccd   = Seq(seq_s, IUPAC.unambiguous_dna)

def LongestSeq(seqlist):
    max_len = 0
    max_seq = None
    for seq in seqlist:
        if len(seq) > max_len:
            max_len = len(seq)
            max_seq = seq
    return max_seq

class ProtCDS:
    def __init__(self):
        self.CDS_FILE = "../Data/CCDS_nucleotide.20150512.fna"
        self.UNIPROT_F= "uniprot_protein.txt"
        self.cds_dict = dict()
        self.uid_list = []
        self.parse_cds_file()
        self.parse_uniprot_file()

    def parse_cds_file(self):
        for cur_record in SeqIO.parse(self.CDS_FILE, "fasta"):
            u_id = cur_record.id.split("|")[0]
            seq  = Seq(cur_record.seq.tostring(), IUPAC.unambiguous_dna)
            self.cds_dict[u_id] = seq

    def getcds(self, cds_id):
        try:
            return self.cds_dict[cds_id]
        except:
            return None

    def parse_uniprot_file(self):
        for line in open(self.UNIPROT_F):
            u_id = line.strip()
            if u_id[0].isalpha():
                self.uid_list.append(u_id)

    def protein_is_in_list(self, uid):
        if uid in self.uid_list:
            return True
        else:
            return False


def PrintCCDSPair(ccd_pair):
    fileout = "tmp"
    fileobj = open(fileout, "w")
    for pair in ccd_pair:
        line = "\t".join(pair)
        fileobj.write(line + "\n")
    fileobj.close()


def PrintMutationPair(ccd_pair):
    fileout = "mut_tmp"
    fileobj = open(fileout, "w")
    for pair in ccd_pair:
        uid  = pair[0]
        seq  = pair[1]
        idx  = 1
        #if uid == "Q9HAN9":
        #    print seq
        for nu3 in [seq[i:(i+3)] for i in xrange(0, len(seq), 3)]:
            #print nu3
            orig_aa, mutate_aa = GenerateAllPossibleMutation(nu3)
            if orig_aa == "*":
                continue
            for each_mutate in mutate_aa:
                content = [uid, idx, nu3, seq3(orig_aa).upper(), seq3(each_mutate).upper()]
                line = "\t".join(map(str, content))
                fileobj.write(line + "\n")
            idx = idx + 1
    fileobj.close()


def GenerateAllPossibleMutation(nucleo_3):
    nucleotide = ["A", "T", "C", "G"]
    mutate_aa  = []
    orig_aa    = Seq("".join(nucleo_3), IUPAC.unambiguous_dna).translate().tostring()
    for i in range(len(nucleo_3)):
        for j in range(len(nucleotide)):
            new3 = list(nucleo_3)
            new3[i] = nucleotide[j]
            aa = Seq("".join(new3), IUPAC.unambiguous_dna).translate().tostring()
            if not aa in mutate_aa and aa  != orig_aa and aa != "*":
                mutate_aa.append(aa)
    return orig_aa, mutate_aa

def main():
    #uniprot2others = "../Data/HUMAN_9606_idmapping.dat"
    ccd2uniprot = "../Data/CCDS2UniProtKB.20150512.txt"
    old_uniprot = ""
    cur_seqlist = []   # current CCD sequence list
    prot_ccd_pair  = []
    all_prot    = []   # keep track of redundant proteins
    cds = ProtCDS()
    ig_flag = 1
    for line in open(ccd2uniprot):
        if ig_flag:
            ig_flag = 0
            continue
        ccd_id, t_type, u_id = line.strip().split("\t")
        if "-" in u_id:
            u_id = u_id.split("-")[0]

        #if u_id != "Q9HAN9":
        #    continue

        if not cds.protein_is_in_list(u_id):
            continue

        if u_id != old_uniprot:
            # start to deal with a new UniProt ID
            if cur_seqlist:
                ccd_seq = LongestSeq(cur_seqlist)
                if not u_id in all_prot:
                    prot_ccd_pair.append((old_uniprot, ccd_seq))
                    all_prot.append(u_id)
            cur_seqlist = []
            ccd_seq     = cds.getcds(ccd_id)
            old_uniprot = u_id
            if ccd_seq is None:
                continue
            else:
                cur_seqlist.append(ccd_seq.tostring())
        else:
            ccd_seq     = cds.getcds(ccd_id)
            if ccd_seq is None:
                continue
            else:
                cur_seqlist.append(ccd_seq.tostring())
    #PrintCCDSPair(prot_ccd_pair)
    PrintMutationPair(prot_ccd_pair)


def test():
    cds = ProtCDS()
    print cds.getcds("CCDS13339.1")

def test2():
    GenerateAllPossibleMutation("AAA")

if __name__ == "__main__":
    main()
    #test2()
