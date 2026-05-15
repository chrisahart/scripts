import numpy as np
import matplotlib.pyplot as plt

homo_s1 = np.array([0.22512488,      0.22940598,      0.23006927,      0.23733470, 0.32162343])
homo_s2 = np.array([0.22505900,      0.22938761,      0.23004975,      0.23729886])
lumo_s1 = np.array([0.34915574])
lumo_s2 = np.array([0.34917847])

ref = np.max(homo_s2)
homo_s1 = homo_s1 - ref
homo_s2 = homo_s2 - ref
lumo_s1 = lumo_s1 - ref
lumo_s2 = lumo_s2 - ref

plt.figure(figsize=(8, 5))
for i, y in enumerate(homo_s1):
    plt.hlines(y, 0.1, 0.9, colors="tab:blue", linewidth=2, label="HOMO S1" if i == 0 else None)
for i, y in enumerate(lumo_s1):
    plt.hlines(y, 0.1, 0.9, colors="tab:green", linewidth=2, label="LUMO S1" if i == 0 else None)
for i, y in enumerate(homo_s2):
    plt.hlines(y, 1.1, 1.9, colors="tab:orange", linewidth=2, label="HOMO S2" if i == 0 else None)
for i, y in enumerate(lumo_s2):
    plt.hlines(y, 1.1, 1.9, colors="tab:red", linewidth=2, label="LUMO S2" if i == 0 else None)

plt.xticks([0.5, 1.5], ["S1", "S2"])
plt.xlabel("State")
plt.ylabel("Energy (shifted, max HOMO S2 = 0)")
plt.title("HOMO/LUMO Energies")
plt.grid(True, linestyle="--", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()
