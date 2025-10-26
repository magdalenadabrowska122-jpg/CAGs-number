def analiza_sekwencji(plik, motyw, ilosc=None):
    import os
    import pandas as pd
    import re
    import numpy as np
    import matplotlib.pyplot as plt

    filename = plik
    with open(filename, mode='r', encoding='utf-8') as f:
        sequence = f.read().upper()
    sequence = sequence.replace("\n", "")

    count = sequence.count(motyw)
    print(f"\nMotyw {motyw} występuje {count} razy w sekwencji.")

    # likalizacja motywu w prostszy sposób,a poniżej z użyciem biblioteki:
    # pozycje = []
    # start = 0
    #
    # while True:
    #     idx = sekwencja.find(motyw, start)
    #     if idx == -1:
    #         break
    #     pozycje.append(idx + 1)
    #     start = idx + 1
    #
    # print(f"Motyw {motyw} występuje {len(pozycje)} razy.")
    # print("Pierwsze pozycje:", pozycje[:10])

    for m in re.finditer(motyw, sequence.upper()):
        start = m.start() + 1
        end = m.end()
        print(f"Motyw {motyw} na pozycji {start}-{end}")

    # Podział na segmenty w prostszy sposób i przy użyciu biblioteki poniżej:
    # k = 30
    # segments = []
    # for i in range(0, len(sequence), k):
    #     fragment = sequence[i:i+k]
    #     segments.append(fragment)
    #
    #
    # for numer, fragment in enumerate(segments, start=1):
    #     ile_razy = fragment.count(motyw)
    #     print(f"Fragment {numer} ({len(fragment)} liter): {ile_razy} razy '{motyw}'")

    segments = []
    k = 30
    for i in range(0, len(sequence), k):
        fragment = sequence[i:i+k]   # wytnij fragment od i do i+k
        segments.append(fragment)
    segments_array = np.array(segments)
    print("Segmenty sekwencji:", segments_array)

    # Znajdź wszystkie pozycje wystąpień motywu
    positions = [m.start() + 1 for m in re.finditer(motyw, sequence.upper())]

    # Wykres słupkowy
    plt.figure(figsize=(100, 2))
    plt.bar(positions, [1]*len(positions), width=2, color='skyblue', edgecolor='black')

    # Opisy osi i tytuł
    plt.xlabel("Pozycja w sekwencji (nt)")
    plt.ylabel("Motyw obecny w sekwencji")
    plt.title(f"Rozmieszczenie motywu '{motyw}' w sekwencji")
    plt.yticks([])  # ukrywamy oś Y, bo nie chcemy jej nazywać. Nic nie wniesie.
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.show()

    df = pd.DataFrame({
        "Position": positions,
        "Segment": [p // k + 1 for p in positions]
    })

    # zapis do CSV
    csv_filename = os.path.join(os.getcwd(), "wyniki.csv")
    df.to_csv(csv_filename, index=False, sep=';')  # średnik dla Excela żeby rozdzielić to na dwie kolumny
    print(f"Plik zapisany! Sprawdź: {os.path.abspath(csv_filename)}")

import tkinter as tk
from tkinter import filedialog, messagebox

def wybierz_sekwencje():
    plik = filedialog.askopenfilename(
        title="Wybierz plik z sekwencją DNA",
        filetypes=[("Pliki tekstowe", "*.txt"), ("Pliki FASTA", "*.fasta"), ("Wszystkie pliki", "*.*")]
    )
    sciezka.set(plik)

def uruchom_analize():
    wybrany_plik = sciezka.get()
    motyw_dna = pole_motyw.get().upper().strip()
    ile_wystapien = pole_ilosc.get().strip()

    if not wybrany_plik or not motyw_dna:
        messagebox.showwarning("Brak danych", "Proszę wybrać plik i wpisać motyw DNA!")
        return

    # liczba wystapien, zamiana na int
    ilosc = int(ile_wystapien) if ile_wystapien.isdigit() else None
    analiza_sekwencji(wybrany_plik, motyw_dna, ilosc)


# OKNO GUI
okno = tk.Tk()
okno.title("Analiza motywów DNA")

sciezka = tk.StringVar()

tk.Button(okno, text="Wybierz sekwencję", command=wybierz_sekwencje).pack(pady=8)
tk.Label(okno, text="Motyw DNA do analizy:").pack()
pole_motyw = tk.Entry(okno)
pole_motyw.pack(pady=5)

tk.Label(okno, text="Ile pierwszych wystąpień motywu analizować (puste = wszystkie):").pack()
pole_ilosc = tk.Entry(okno)
pole_ilosc.pack(pady=5)

tk.Button(okno, text="Rozpocznij analizę", command=uruchom_analize).pack(pady=8)

okno.mainloop()