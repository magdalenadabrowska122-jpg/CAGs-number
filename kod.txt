import os
import pandas as pd
filename = ('sekwencja.txt')
with open(filename, mode = 'r', encoding = 'utf-8') as f:
    sequence = f.read().upper()
sequence = sequence.replace("\n", "")
motif = "CAG"
count = sequence.count(motif)
print(f"\nMotyw {motif} występuje {count} razy w sekwencji.")
import re
for m in re.finditer(motif, sequence.upper()):
    start = m.start() + 1
    end = m.end()
    print(f"Motyw {motif} na pozycji {start}-{end}")
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

# Tworzymy wykres słupkowy
plt.figure(figsize=(100, 2))
plt.bar(positions, [1]*len(positions), width=2, color='skyblue', edgecolor='black')

# Opisy osi i tytuł
plt.xlabel("Pozycja w sekwencji (nt)")
plt.ylabel("Motyw obecny")
plt.title(f"Rozmieszczenie motywu '{motif}' w sekwencji")
plt.yticks([])  # ukrywamy oś Y, bo nie niesie informacji
plt.grid(axis='x', linestyle='--', alpha=0.5)
plt.show()
segment_size = k
df = pd.DataFrame({
    "Position": positions,
    "Segment": [p // segment_size + 1 for p in positions]
})

# zapis do CSV
csv_filename = os.path.join(os.getcwd(), "wyniki.csv")
df.to_csv(csv_filename, index=False, sep=';')  # średnik dla Excela
print(f"Plik zapisany! Sprawdź: {os.path.abspath(csv_filename)}")