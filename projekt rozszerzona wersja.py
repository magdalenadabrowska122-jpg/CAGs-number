from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import os
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import tkinter as tk
from tkinter import filedialog, messagebox

Entrez.email = "magdalenadabrowska122@gmail.com"


#Analiza sekwencji z NCBI
def analiza_sekwencji(genbank_id, motyw):
    # pobranie sekwencji
    a = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
    b = SeqIO.read(a, "fasta")
    a.close()

    print("ID rekordu:", b.id)
    print("Opis:", b.description)
    print("Pierwsze 50 nukleotydów:", b.seq[:50])
    print("Długość całej sekwencji:", len(b.seq))

    # Szukanie motywu
    pozycje = []
    start = 0
    while True:
        idx = b.seq.find(motyw, start)
        if idx == -1:
            break
        pozycje.append(idx)
        start = idx + 1

    print(
        f"\nMotyw {motyw} występuje {len(pozycje)} razy na pozycjach: {pozycje[:10]}{'...' if len(pozycje) > 10 else ''}")
# Zabezpieczenie przed brakiem potencjlanego wyniku
    if not pozycje:
        print(f"Motyw {motyw} nie występuje w sekwencji!")
        messagebox.showinfo("Brak wyników", f"Motyw {motyw} nie został znaleziony w sekwencji.")
        return

    # Wykres matplotlib
    fig, ax = plt.subplots(figsize=(12, 2))
    ax.plot([pos for pos in pozycje], [0] * len(pozycje), 'ro', markersize=5)
    ax.set_title(f"Motyw '{motyw}' w sekwencji {b.id}")
    ax.set_xlabel("Pozycja nukleotydu")
    ax.set_yticks([])
    ax.grid(axis='x', linestyle='--', alpha=0.7)

    pdf_filename = os.path.join(os.getcwd(), f"wykres_{b.id}_{motyw}.pdf")
    with PdfPages(pdf_filename) as pdf:
        pdf.savefig(fig)
    print(f"Plik PDF z wykresem zapisany! Sprawdź: {pdf_filename}")
    plt.show()

    csv_filename1 = os.path.join(os.getcwd(), f"pozycje_{b.id}_{motyw}.csv")
    pd.DataFrame({"Pozycja": pozycje}).to_csv(csv_filename1, index=False, sep=';')
    print(f"Plik CSV zapisany! Sprawdź: {csv_filename1}")

    # Wykres interaktywny
    fig = go.Figure()
    for pos in pozycje:
        fig.add_trace(go.Scatter(
            x=[pos, pos + len(motyw)],
            y=[0, 0],
            mode="lines",
            line=dict(color="red", width=10),
            name=motyw,
            hovertemplate=f"Motyw: {motyw}<br>Pozycja: {pos}"
        ))
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(
        title=f"Motyw '{motyw}' w sekwencji {b.id}",
        xaxis_title="Pozycja nukleotydu",
        showlegend=True,
        height=200
    )
    fig.show()


#GUI
def uruchom_analize():
    genbank_id = pole_genbank.get().strip()
    motyw = pole_motyw.get().upper().strip()

    if not genbank_id or not motyw:
        messagebox.showwarning("Brak danych", "Proszę wpisać GenBank ID i motyw DNA!")
        return

    analiza_sekwencji(genbank_id, motyw)


okno = tk.Tk()
okno.title("Analiza motywów DNA z NCBI")

tk.Label(okno, text="GenBank ID:").pack()
pole_genbank = tk.Entry(okno)
pole_genbank.pack(pady=5)

tk.Label(okno, text="Motyw DNA:").pack()
pole_motyw = tk.Entry(okno)
pole_motyw.pack(pady=5)

tk.Button(okno, text="Rozpocznij analizę", command=uruchom_analize).pack(pady=10)

okno.mainloop()