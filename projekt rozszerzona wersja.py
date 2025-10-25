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
# OK. Tworzenie wykresu wystepowania danych motywów
import matplotlib.pyplot as plt
kolory= {"CAG": "red", "ATG": "yellow", "CCC": "violet"}
plt.figure(figsize=(12,2))
plt.title("Motywy w sekwencji TP53")
plt.xlabel("pozyccja nukleotydu")
plt.yticks([])
for motyw, pozycje in wyniki.items():
    for pos in pozycje:
        plt.plot([pos, pos+len(motyw)],[0, 0], color=kolory[motyw], linewidth=5)

handles = [plt.Line2D([0],[0], color=kolory[m], lw=5) for m in motywy]
plt.legend(handles, motywy)
plt.show()
print("wygenerowany wykres")
import plotly.graph_objects as go

motywy = ["CAG", "ATG", "CCC"]
kolory = {"CAG": "red", "ATG": "green", "CCC": "blue"}

fig = go.Figure()

for motyw in motywy:
    pozycje = wyniki[motyw]
    for pos in pozycje:
        fig.add_trace(go.Scatter(
            x=[pos, pos+len(motyw)],
            y=[0, 0],
            mode="lines",
            line=dict(color=kolory[motyw], width=10),
            name=motyw,
            hovertemplate=f"Motyw: {motyw}<br>Pozycja: {pos}"
        ))


fig.update_yaxes(showticklabels=False)
fig.update_layout(
    title="Motywy w sekwencji TP53",
    xaxis_title="Pozycja nukleotydu",
    showlegend=True,
    height=200
)

fig.show()
print("wygenerowany wykres interaktywny")
