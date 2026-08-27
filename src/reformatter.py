# Reformatter for txt refractive index data
# Use when it is not organized in 3 columns
# TO DO: Add the ability to choose the file to reformat by the user

from pathlib import Path
import pandas as pd

src = Path("/mnt/data/Au.txt")
text = src.read_text()

lines = [l.strip() for l in text.splitlines()]

n_data = []
k_data = []
mode = None

for line in lines:
    if not line:
        continue
    if line.lower() == "wl\tn" or line.lower().replace(" ", "") == "wl\tn":
        mode = "n"
        continue
    if line.lower() == "wl\tk" or line.lower().replace(" ", "") == "wl\tk":
        mode = "k"
        continue

    parts = line.split()
    if len(parts) == 2:
        wl, val = map(float, parts)
        if mode == "n":
            n_data.append((wl, val))
        elif mode == "k":
            k_data.append((wl, val))

df_n = pd.DataFrame(n_data, columns=["wl", "n"])
df_k = pd.DataFrame(k_data, columns=["wl", "k"])

df = pd.merge(df_n, df_k, on="wl", how="outer")
out = Path("/mnt/data/Au_reformatted.txt")
df.to_csv(out, sep="\t", index=False)

print(str(out))
