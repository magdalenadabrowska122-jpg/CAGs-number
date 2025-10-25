from Bio import Entrez, SeqIO


Entrez.email = "magdalenadabrowska122@gmail.com"
# z NCBI pobrałam id geny TP53
genbank_id = "NM_000546"
a = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
b = SeqIO.read(a, "fasta")
a.close()
print("ID rekordu:", b.id)
print("Opis:", b.description)
print("Pierwsze 50 nukleotydów:", b.seq[:50])
print("Długość całej sekwencji:", len(b.seq))
# Wyszukanie wielu motywów jednoczenie. Wybrałam motywy CAG, CCC, ATG
motywy = ["CAG", "ATG", "CCC"]
# Postaram się użyć słownika
wyniki = {}
for motyw in motywy:
    pozycje = []
    start = 0
    while True:
        idx = b.seq.find(motyw, start)
        if idx == -1:
            break
        pozycje.append(idx)
        start = idx + 1
    wyniki[motyw] = pozycje

for motyw, pozycje in wyniki.items():
    print(f"Motyw {motyw} występuje {len(pozycje)} razy na pozycjach: {pozycje[:10]}{'...' if len(pozycje) > 10 else ''}")