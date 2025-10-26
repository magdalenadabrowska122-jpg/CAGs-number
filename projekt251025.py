import os
import pandas as pd
filename = ('sekwencja.txt')
with open(filename, mode = 'r', encoding = 'utf-8') as f:
    sequence = f.read().upper()
sequence = sequence.replace("\n", "")
motif = "CAG"
count = sequence.count(motif)
print(f"\nMotyw {motif} występuje {count} razy w sekwencji.")
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
import re
for m in re.finditer(motif, sequence.upper()):
    start = m.start() + 1
    end = m.end()
    print(f"Motyw {motif} na pozycji {start}-{end}")
# Podział na segmenty w prostszy sposób i przy użyciu biblioteki poniżej:
# k = 30
# segments = []
# for i in range(0, len(sequence), k):
#     fragment = sequence[i:i+k]
#     segments.append(fragment)
#
#
# for numer, fragment in enumerate(segments, start=1):
#     ile_razy = fragment.count(motif)
#     print(f"Fragment {numer} ({len(fragment)} liter): {ile_razy} razy '{motif}'")

import numpy as np
segments = []
k = 30
for i in range(0, len(sequence), k):
    fragment = sequence[i:i+k]   # wytnij fragment od i do i+k
    segments.append(fragment)
segments_array = np.array(segments)
print("Segmenty sekwencji:", segments_array)
import re
import matplotlib.pyplot as plt

# Znajdź wszystkie pozycje wystąpień motywu
positions = [m.start() + 1 for m in re.finditer(motif, sequence.upper())]

# Wykres słupkowy
plt.figure(figsize=(100, 2))
plt.bar(positions, [1]*len(positions), width=2, color='skyblue', edgecolor='black')

# Opisy osi i tytuł
plt.xlabel("Pozycja w sekwencji (nt)")
plt.ylabel("Motyw obecny w sekwencji")
plt.title(f"Rozmieszczenie motywu '{motif}' w sekwencji")
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
