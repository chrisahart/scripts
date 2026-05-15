import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.transformations import set_dimensions
from MDAnalysis.lib.distances import distance_array, apply_PBC
import warnings
warnings.filterwarnings('ignore')

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': '14',
          'axes.titlesize': '16',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

# ─────────────────────────────────────────────
# Generated with help of Claude
# ─────────────────────────────────────────────

# folder = '/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/100/md/zeng_input/optb88'
# BOX        = [13.9553247, 8.90701634, 31.33880111, 90.0, 90.0, 90.0]

# folder = '/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/110/md/pbe-d-d3-300k-temptol-30-nose-100'
# BOX        = [13.155592560000001, 8.9070163400000002, 31.389484599999999, 90.0, 90.0, 90.0]
# TRAJ       = '{}/tio2-pos-1.xyz'.format(folder)
# TRAJ       = '{}/tio2-pos-1_center.xyz'.format(folder)
# TRAJ_START = 0

# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/archer/ahart/110/water/md/pbe-frozen-tio2/sep-h-3.5'
# BOX        = [12.9824805,  17.76,       33.3469092, 90.0, 90.0, 90.0]
# TRAJ       = '{}/system.xyz'.format(folder)
# TRAJ       = '{}/tio2-pos-1.xyz'.format(folder)
# TRAJ_START = 0

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/archer/ahart/110/water/md/pbe-frozen-tio2/sep-o-3.2-center'
BOX        = [12.9824805,  17.76,       32.239747139999999, 90.0, 90.0, 90.0]
folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/archer/ahart/110/water/md/pbe-frozen-tio2/sep-o-3.4-center'
BOX        = [12.9824805,  17.76,       32.639747139999997, 90.0, 90.0, 90.0]
folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/archer/ahart/110/water/md/pbe-frozen-tio2/sep-o-3.6-center'
BOX        = [12.9824805,  17.76,       33.039747140000003, 90.0, 90.0, 90.0]
TRAJ       = '{}/system.xyz'.format(folder)
TRAJ_START = 0
BULK_Z_MIN = 17
BULK_Z_MAX = 27
TRAJ       = '{}/tio2-pos-1.xyz'.format(folder)
TRAJ_START = 700
BULK_Z_MIN = 17
BULK_Z_MAX = 27

STRIDE = 1
save_fig = True
OUTDIR = folder

# z is the surface-normal direction (longest box vector)
SURF_NORMAL = 2  # index: 0=x, 1=y, 2=z
box3 = np.array(BOX[:3])
Lz = BOX[SURF_NORMAL]

# Density profile binning
BIN_WIDTH = 0.05  # Å — bin width along z
N_BINS = int(Lz / BIN_WIDTH)
BIN_EDGES = np.linspace(0, Lz, N_BINS + 1)
BIN_CENTERS = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])

# O-H bond cutoff for identifying water O vs TiO2 O
OH_CUTOFF = 1.30  # Å — generous to catch all water O

print("=" * 60)
print("TiO2/Water Interface MD Analysis")
print(f"Output folder  : {OUTDIR}")
print(f"Box            : {BOX[:3]}")
print(f"TRAJ_START     : frame {TRAJ_START}")
print("=" * 60)

# ─────────────────────────────────────────────
# LOAD
# ─────────────────────────────────────────────
u = mda.Universe(TRAJ, format='XYZ')
u.trajectory.add_transformations(set_dimensions(BOX))

n_frames_total = len(u.trajectory)
n_frames = (n_frames_total - TRAJ_START) // STRIDE
prod_slice = slice(TRAJ_START, n_frames_total, STRIDE)

print(f"Total frames   : {n_frames_total}")
print(f"Skipped        : {TRAJ_START} frames (equilibration)")
print(f"Stride         : {STRIDE}")
print(f"Frames used    : {n_frames}")

# Atom selections
# In CP2K XYZ: Ti = titanium, O = all oxygens, H = hydrogens
Ti_sel = u.select_atoms('name Ti')
O_sel = u.select_atoms('name O')
H_sel = u.select_atoms('name H')

print(f"\nAtom counts:")
print(f"  Ti : {len(Ti_sel)}")
print(f"  O  : {len(O_sel)}  (TiO2 + water)")
print(f"  H  : {len(H_sel)}  (water only)")
print(f"  Water molecules (expected) : {len(H_sel) // 2}")

# ─────────────────────────────────────────────
# IDENTIFY WATER O vs TiO2 O
# Strategy: water O has at least one H within OH_CUTOFF
# Do this once on first production frame — atom types are fixed
# ─────────────────────────────────────────────
print("\nIdentifying water O vs TiO2 O...")

u.trajectory[TRAJ_START]
dims = u.trajectory.ts.dimensions
pos_O = O_sel.positions.copy()
pos_H = H_sel.positions.copy()

# distance_array from MDAnalysis uses PBC correctly
dist_OH = distance_array(pos_O, pos_H, box=dims)  # (N_O, N_H)
is_water_O = np.any(dist_OH < OH_CUTOFF, axis=1)  # True if bonded to any H

water_O_indices = O_sel.indices[is_water_O]
tio2_O_indices = O_sel.indices[~is_water_O]

waterO_sel = u.select_atoms('index ' + ' '.join(map(str, water_O_indices)))
tio2O_sel = u.select_atoms('index ' + ' '.join(map(str, tio2_O_indices)))

print(f"  Water O atoms  : {len(waterO_sel)}")
print(f"  TiO2 O atoms   : {len(tio2O_sel)}")
print(f"  Water molecules: {len(waterO_sel)}")

# ─────────────────────────────────────────────
# FIND Ti LAYER POSITIONS
# ─────────────────────────────────────────────
print("\nLocating Ti layer positions...")

# Collect mean z per Ti atom over production frames
ti_z_sum = np.zeros(len(Ti_sel))
for ts in u.trajectory[prod_slice]:
    ti_z_sum += Ti_sel.positions[:, SURF_NORMAL]
ti_z_mean_pos = ti_z_sum / n_frames  # mean z per Ti atom

# Find discrete Ti layers by greedy clustering:
# Start from the highest Ti z, assign all Ti within LAYER_TOL to the same
# layer and compute their mean z, then repeat on the remaining atoms.
LAYER_TOL = 1.0  # Å
remaining = ti_z_mean_pos.copy()
layers = []
while len(remaining) > 0:
    top = remaining.max()
    mask = np.abs(remaining - top) < LAYER_TOL
    layers.append(float(remaining[mask].mean()))
    remaining = remaining[~mask]
layers = np.sort(layers)  # ascending z

print(f"  Ti layers found : {len(layers)}")
for i, lz in enumerate(layers):
    print(f"    Layer {i + 1}: z = {lz:.4f} Å")

# Outermost faces: bottom = min layer z, top = max layer z
z_ref_bottom = layers[0]  # faces water below (smaller z)
z_ref_top = layers[-1]  # faces water above (larger z)

print(f"  h=0 bottom face : {z_ref_bottom:.4f} Å")
print(f"  h=0 top face    : {z_ref_top:.4f} Å")
print(f"  Slab thickness  : {z_ref_top - z_ref_bottom:.4f} Å")

# ─────────────────────────────────────────────
# WATER DENSITY PROFILE
# Bin water O atoms along z, then shift to h from each Ti face
# ─────────────────────────────────────────────
print("\nCalculating water density profile (streaming)...")

hist_raw = np.zeros(N_BINS)
n_counted = 0

waterO_idx = waterO_sel.indices

for ts in u.trajectory[prod_slice]:
    pos = u.atoms.positions.copy()
    pos = apply_PBC(pos, ts.dimensions)
    wO_z = pos[waterO_idx, SURF_NORMAL]
    hist_raw += np.histogram(wO_z, bins=BIN_EDGES)[0]
    n_counted += 1

# Normalise: number density (Å⁻³)
bin_vol = BOX[0] * BOX[1] * BIN_WIDTH  # Å³
rho_bins = hist_raw / (n_counted * bin_vol)

# ─────────────────────────────────────────────
# BUILD h PROFILES FROM EACH SURFACE
# Slab is in the MIDDLE of the box:
#   bottom water: z < z_ref_bottom  → h = z_ref_bottom - z  (h>0 into water)
#   top water:    z > z_ref_top     → h = z - z_ref_top     (h>0 into water)
# ─────────────────────────────────────────────

# Water region extents — account for PBC wrapping.
# Bottom water: z < z_ref_bottom AND z > z_ref_top (wraps across boundary)
# Top water:    z > z_ref_top (main water region)
# Total water thickness on each side estimated from box geometry.
water_bot_thickness = z_ref_bottom + (Lz - z_ref_top)  # wraps PBC
water_top_thickness = Lz - z_ref_top
h_max = min(water_bot_thickness, water_top_thickness) - 0.5

# Bottom surface: water wraps across PBC boundary
# Bins with z < z_ref_bottom contribute h = z_ref_bottom - z  (normal)
# Bins with z > z_ref_top contribute h = z_ref_bottom + (Lz - z) (wrapped)
h_bottom_normal = z_ref_bottom - BIN_CENTERS  # z < z_ref_bottom
h_bottom_wrapped = z_ref_bottom + (Lz - BIN_CENTERS)  # z > z_ref_top (PBC image)

mask_bot_n = (h_bottom_normal >= 0.0) & (h_bottom_normal <= h_max)
mask_bot_w = (h_bottom_wrapped >= 0.0) & (h_bottom_wrapped <= h_max)

h_b = np.concatenate([h_bottom_normal[mask_bot_n], h_bottom_wrapped[mask_bot_w]])
rho_b = np.concatenate([rho_bins[mask_bot_n], rho_bins[mask_bot_w]])
sort_b = np.argsort(h_b)
h_b = h_b[sort_b]
rho_b = rho_b[sort_b]

# Top surface: water has larger z than slab → h = z - z_ref_top
h_top = BIN_CENTERS - z_ref_top
mask_top = (h_top >= -5.0) & (h_top <= h_max)
h_t = h_top[mask_top]
rho_t = rho_bins[mask_top]
sort_t = np.argsort(h_t)
h_t = h_t[sort_t]
rho_t = rho_t[sort_t]

bulk_mask = (BIN_CENTERS >= BULK_Z_MIN) & (BIN_CENTERS <= BULK_Z_MAX)
rho_bulk = rho_bins[bulk_mask].mean() if bulk_mask.sum() > 0 else None

print(f"  Bulk z range              : {BULK_Z_MIN} – {BULK_Z_MAX} Å  ({bulk_mask.sum()} bins)")
print(f"  Bulk water number density : {rho_bulk:.4f} Å\u207b\u00b3" if rho_bulk else "  Not enough bins in bulk range")

# Expected bulk water density ~0.0334 Å⁻³ (1 g/cm³)
rho_ref = 0.03334  # Å⁻³

# ─────────────────────────────────────────────
# FIGURE 1: BOTH SURFACES OVERLAPPED
# h=0 at outermost Ti layer on each face
# ─────────────────────────────────────────────
fig_dens_both, ax_dens_both = plt.subplots(figsize=(7, 5))
ax_dens_both.plot(h_b, rho_b, 'k-', lw=1.5, label='Bottom surface')
ax_dens_both.plot(h_t, rho_t, 'k--', lw=1.5, label='Top surface')
if rho_bulk:
    ax_dens_both.axhline(rho_bulk, color='gray', ls='--', lw=1.2,
                         label=f'Bulk = {rho_bulk:.4f} Å\u207b\u00b3')
ax_dens_both.axhline(rho_ref, color='gray', ls=':', lw=1.2,
                     label=f'Ref = {rho_ref:.4f} Å\u207b\u00b3')
ax_dens_both.axvline(0.0, color='lightgray', ls='-', lw=1, zorder=0)
ax_dens_both.set_xlabel(r'h (Å)  [h=0 at outermost Ti layer]')
ax_dens_both.set_ylabel(r'$\rho_{\mathrm{O}}$ (Å$^{-3}$)')
ax_dens_both.set_xlim([-5, h_max])
ax_dens_both.set_ylim(bottom=0)
# ax_dens_both.legend()
fig_dens_both.tight_layout()
if save_fig: fig_dens_both.savefig('{}/water_density_profile.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# FIGURE 2: FULL z PROFILE (absolute, not shifted)
# shows slab + both water layers
# ─────────────────────────────────────────────
fig_dens_full, ax_dens_full = plt.subplots(figsize=(8, 4))
ax_dens_full.plot(BIN_CENTERS, rho_bins, 'k-', lw=1.2)
if rho_bulk:
    ax_dens_full.axhline(rho_bulk, color='gray', ls='--', lw=1.2,
                         label=f'Bulk = {rho_bulk:.4f} Å\u207b\u00b3')
ax_dens_full.axhline(rho_ref, color='gray', ls=':', lw=1.2,
                     label=f'Ref = {rho_ref:.4f} Å\u207b\u00b3')
ax_dens_full.axvline(z_ref_bottom, color='steelblue', ls='--', lw=1.2,
                     label=f'Ti bottom face = {z_ref_bottom:.2f} Å')
ax_dens_full.axvline(z_ref_top, color='firebrick', ls='--', lw=1.2,
                     label=f'Ti top face = {z_ref_top:.2f} Å')
ax_dens_full.set_xlabel(r'z (Å)')
ax_dens_full.set_ylabel(r'$\rho_{\mathrm{O}}$ (Å$^{-3}$)')
ax_dens_full.set_xlim([0, Lz])
ax_dens_full.set_ylim(bottom=0)
# ax_dens_full.legend(fontsize=11)
fig_dens_full.tight_layout()
if save_fig: fig_dens_full.savefig('{}/water_density_profile_full.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Total frames        : {n_frames_total}")
print(f"Production frames   : {n_frames}  (stride={STRIDE})")
print(f"Water O atoms       : {len(waterO_sel)}")
print(f"TiO2 O atoms        : {len(tio2O_sel)}")
print(f"Ti atoms            : {len(Ti_sel)}")
print(f"Ti layers           : {len(layers)}")
print(f"Ti bottom face (h=0): {z_ref_bottom:.4f} Å")
print(f"Ti top face (h=0)   : {z_ref_top:.4f} Å")
print(f"Slab thickness      : {z_ref_top - z_ref_bottom:.4f} Å")
if rho_bulk:
    print(f"Bulk density        : {rho_bulk:.4f} Å⁻³  (ref: {rho_ref:.4f} Å⁻³)")
print("=" * 60)
print("\nFigures saved:")
for f in ['water_density_profile.png', 'water_density_profile_full.png']:
    print(f"  {os.path.join(OUTDIR, f)}")

if __name__ == "__main__":
    print('Finished.')
    plt.show()