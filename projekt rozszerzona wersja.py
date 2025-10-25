from Bio import Entrez, SeqIO


Entrez.email = "magdalenadabrowska122@gmail.com"
# z NCBI pobrałam id geny TP53
genbank_id = "NM_000546"
a = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
b = SeqIO.read(a, "fasta")
a.close()
