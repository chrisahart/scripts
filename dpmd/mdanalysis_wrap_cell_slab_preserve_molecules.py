import MDAnalysis as mda
import MDAnalysis.transformations as trans
import numpy as np

folder = '/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/110/md/pbe-d-d3-300k-temptol-30-nose-100'
cell = [13.15559256, 8.90701634, 31.3894846, 90.0, 90.0, 90.0]
filename_input = "{}/tio2-pos-1.xyz".format(folder)
filename_output = "{}/tio2-pos-1_center_preserve.xyz".format(folder)
filename_ref = "{}/system.xyz".format(folder)

# Load systems
u = mda.Universe(filename_input, format="XYZ")
ref = mda.Universe(filename_ref, format="XYZ")

u.dimensions = cell
ref.dimensions = cell

# ─────────────────────────────────────────────
# Ti-only alignment groups
# ─────────────────────────────────────────────
mobile = u.select_atoms("name Ti")
reference = ref.select_atoms("name Ti")

if len(mobile) != len(reference):
    raise ValueError("Mismatch in Ti atoms between trajectory and reference")

# apply alignment transform
u.trajectory.add_transformations(
    trans.fit_rot_trans(mobile, reference)
)

# ─────────────────────────────────────────────
# WATER INDEXING (H, H, O ordering)
# ─────────────────────────────────────────────
start = 144  # first water atom index (0-based)
n_atoms = len(u.atoms)

if (n_atoms - start) % 3 != 0:
    raise ValueError("Water region not divisible by 3 — check indexing")

water_idx = np.arange(start, n_atoms).reshape(-1, 3)  # (n_mol, 3)

# ─────────────────────────────────────────────
# TRAJECTORY LOOP
# ─────────────────────────────────────────────
with open(filename_output, "w") as f:

    for ts in u.trajectory:

        coords = u.atoms.positions.copy()
        box = u.dimensions[:3]

        # --- molecule-aware wrapping (water only) ---
        mol_coords = coords[water_idx]        # (n_mol, 3, 3)

        # center of geometry (robust for HHO)
        com = mol_coords.mean(axis=1)         # (n_mol, 3)

        # wrap COM into box
        com_wrapped = com - box * np.round(com / box)

        shift = com_wrapped - com             # (n_mol, 3)

        # apply translation to whole molecule
        mol_coords += shift[:, None, :]

        # write back
        coords[water_idx] = mol_coords

        # --- optional: wrap everything into primary box ---
        coords = coords - box * np.floor(coords / box)

        # --- WRITE FRAME ---
        f.write(f"{len(coords)}\n\n")
        for atom, coord in zip(u.atoms.names, coords):
            f.write(f"{atom} {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f}\n")

print(f"Finished. Output written to: {filename_output}")