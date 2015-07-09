from Bio import Entrez

rec = Entrez.read(Entrez.esearch(db="nucleotide", term="NM_003404.4"))

p_handle = Entrez.efetch(db="protein", id=rec["IdList"][0], rettype="fasta")

print p_handle.read()

from Bio import SeqIO
record = SeqIO.read(open("single.fasta"), "fasta")
