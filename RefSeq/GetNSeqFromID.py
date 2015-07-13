'''
  Get the nucleotide sequence for each protein
'''

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

Entrez.email = 'ajingnk@gmail.com'



def GetSeqByProtID(protid):
    rec = Entrez.read(Entrez.esearch(db="protein", term= protid))
    p_handle = Entrez.efetch(db="protein", id=rec["IdList"][0], rettype="fasta")
    for cur_record in SeqIO.parse(p_handle, "fasta"):
        u_id = cur_record.id
        m_dna= cur_record.seq.tostring()
    return u_id, m_dna

def GetSeqByRefID(refid):
    rec = Entrez.read(Entrez.esearch(db="nucleotide", term= refid))
    p_handle = Entrez.efetch(db="nucleotide", id=rec["IdList"][0], rettype="fasta")
    for cur_record in SeqIO.parse(p_handle, "fasta"):
        u_id = cur_record.id
        m_dna= Seq(cur_record.seq.tostring(), IUPAC.unambiguous_dna)
        #print m_dna.transcribe().translate()
        #print dir(m_dna.transcribe())
        #print m_dna.translate()
    return u_id, m_dna


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

def GetCDS(u_id):


def main():
    uniprot2refseq = "../Data/HUMAN_9606_idmapping.dat"
    cur_uniprot = ""
    cur_seqlist = []   # current CCD sequence list
    prot_ccd_pair  = []
    for line in open(uniprot2refseq):
        u_id, t_type, t_id = line.strip().split("\t")
        if t_type == "CCDS":
            if u_id != cur_uniprot:
                # start to deal with a new UniProt ID
                if cur_seqlist:
                    ccd_seq = LongestSeq(cur_seqlist)
                    prot_ccd_pair.append((u_id, ccd_seq))
                cur_seqlist = []
                ccd_seq     = GetCDS(u_id)
                cur_seqlist.append(ccd_seq)
            else:
                ccd_seq     = GetCDS(u_id)
                cur_seqlist.append(ccd_seq)


def test():
    print GetSeqByRefID("NM_006761.4")
    print GetSeqByProtID("NP_006752.1")


if __name__ == "__main__":
    test()
    seq = """
    ATGACAATGGATAAAAGTGAGCTGGTACAGAAAGCCAAACTCGCTGAGCAGGCTGAGCGA
    TATGATGATATGGCTGCAGCCATGAAGGCAGTCACAGAACAGGGGCATGAACTCTCCAAC
    GAAGAGAGAAATCTGCTCTCTGTTGCCTACAAGAATGTGGTAGGCGCCCGCCGCTCTTCC
    TGGCGTGTCATCTCCAGCATTGAGCAGAAAACAGAGAGGAATGAGAAGAAGCAGCAGATG
    GGCAAAGAGTACCGTGAGAAGATAGAGGCAGAACTGCAGGACATCTGCAATGATGTTCTG
    GAGCTGTTGGACAAATATCTTATTCCCAATGCTACACAACCAGAAAGTAAGGTGTTCTAC
    TTGAAAATGAAAGGAGATTATTTTAGGTATCTTTCTGAAGTGGCATCTGGAGACAACAAA
    CAAACCACTGTGTCGAACTCCCAGCAGGCTTACCAGGAAGCATTTGAAATTAGTAAGAAA
    GAAATGCAGCCTACACACCCAATTCGTCTTGGTCTGGCACTAAATTTCTCAGTCTTTTAC
    TATGAGATTCTAAACTCTCCTGAAAAGGCCTGTAGCCTGGCAAAAACGGCATTTGATGAA
    GCAATTGCTGAATTGGATACGCTGAATGAAGAGTCTTATAAAGACAGCACTCTGATCATG
    CAGTTACTTAGGGACAATCTCACTCTGTGGACATCGGAAAACCAGGGAGACGAAGGAGAC
    GCTGGGGAGGGAGAGAACTAA
    """
    GetTranslatedSeq(seq)
