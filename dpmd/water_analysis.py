"""
Water Box MD Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────
TRAJ      = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/cubic/64/md/pbe/md/from-tip3p-1-ns/tio2-pos-1.xyz'
BOX = [12.4, 12.4, 12.4, 90.0, 90.0, 90.0]
RDF_RMAX = 12.4/2
N_WATERS = 64

# TRAJ      = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/cubic/128/md/pbe/from-tip3p-1-ns/tio2-pos-1.xyz'
# BOX = [15.65, 15.65, 15.65, 90.0, 90.0, 90.0]
# RDF_RMAX = 15.65/2
# N_WATERS = 128

TIMESTEP = 1e-3  # ps
STRIDE = 1  # read every Nth frame — increase for speed
RDF_BINS = 200
HBOND_DIST = 3.5  # Å  O...O cutoff
HBOND_ANG = 30.0  # degrees

# ─────────────────────────────────────────────
# LOAD WITH MDANALYSIS (streaming, low memory)
# ─────────────────────────────────────────────
print(f"\nLoading: {TRAJ}")

u = mda.Universe(TRAJ, format='XYZ')

# Set box for all frames via a transformation
from MDAnalysis.transformations import set_dimensions

transform = set_dimensions(BOX)
u.trajectory.add_transformations(transform)

n_frames_total = len(u.trajectory)
n_frames = n_frames_total // STRIDE
times = np.arange(n_frames) * TIMESTEP * STRIDE

print(f"Total frames  : {n_frames_total}")
print(f"Stride        : {STRIDE}")
print(f"Frames used   : {n_frames}")
print(f"Total time    : {n_frames * TIMESTEP * STRIDE:.2f} ps")

# Atom selections
O_sel = u.select_atoms('name O')
H_sel = u.select_atoms('name H')
print(f"O atoms       : {len(O_sel)}")
print(f"H atoms       : {len(H_sel)}")

# ─────────────────────────────────────────────
# 1. RDFs via MDAnalysis InterRDF (C-level, fast)
# ─────────────────────────────────────────────
print("\nCalculating RDFs (MDAnalysis InterRDF)...")

rdf_OO = InterRDF(O_sel, O_sel, nbins=RDF_BINS, range=(0.5, RDF_RMAX), exclusion_block=(1, 1))
rdf_OH = InterRDF(O_sel, H_sel, nbins=RDF_BINS, range=(0.5, RDF_RMAX))
rdf_HH = InterRDF(H_sel, H_sel, nbins=RDF_BINS, range=(0.5, RDF_RMAX), exclusion_block=(1, 1))

rdf_OO.run(step=STRIDE, verbose=True)
rdf_OH.run(step=STRIDE, verbose=True)
rdf_HH.run(step=STRIDE, verbose=True)

r_OO = rdf_OO.results.bins
g_OO = rdf_OO.results.rdf
r_OH = rdf_OH.results.bins
g_OH = rdf_OH.results.rdf
r_HH = rdf_HH.results.bins
g_HH = rdf_HH.results.rdf

# First peak and minimum O-O
dr = r_OO[1] - r_OO[0]
peak_idx = int(2.0 / dr) + np.argmax(g_OO[int(2.0 / dr):int(4.0 / dr)])
min_idx = peak_idx + np.argmin(g_OO[peak_idx:peak_idx + int(2.0 / dr)])
print(f"  O-O first peak : {r_OO[peak_idx]:.3f} Å  (ref: 2.80 Å)")
print(f"  O-O first min  : {r_OO[min_idx]:.3f} Å  (ref: 3.50 Å)")

# Coordination number
rho_O = len(O_sel) / (BOX[0] * BOX[1] * BOX[2])
try:
    coord_num = 4 * np.pi * rho_O * np.trapezoid(g_OO[:min_idx] * r_OO[:min_idx] ** 2, r_OO[:min_idx])
except AttributeError:
    coord_num = 4 * np.pi * rho_O * np.trapz(g_OO[:min_idx] * r_OO[:min_idx] ** 2, r_OO[:min_idx])
print(f"  Coordination N : {coord_num:.2f}  (ref: 4.5)")

# ─────────────────────────────────────────────
# 2. MSD — vectorized, streaming
# ─────────────────────────────────────────────
print("\nCalculating MSD (streaming)...")

msd = np.zeros(n_frames)
ref = None
box3 = np.array(BOX[:3])

for fi, ts in enumerate(u.trajectory[::STRIDE]):
    pos_O = O_sel.positions.copy()
    if fi == 0:
        ref = pos_O.copy()
    disp = pos_O - ref
    disp -= np.round(disp / box3) * box3
    msd[fi] = np.mean(np.sum(disp ** 2, axis=1))

# Linear fit middle 60% of trajectory
s = int(0.2 * n_frames)
e = int(0.8 * n_frames)
slope, intercept = np.polyfit(times[s:e], msd[s:e], 1)
D = slope / 6.0
D_SI = D * 1e-20 / 1e-12
print(f"  Diffusion coeff : {D:.4f} Å²/ps = {D_SI:.2e} m²/s")
print(f"  (Exp TIP3P ~5.0e-9, Exp water ~2.3e-9 m²/s)")

# ─────────────────────────────────────────────
# 3. H-BONDS, O-H LENGTHS, H-O-H ANGLES, ENERGY
#    — single streaming pass
# ─────────────────────────────────────────────
print("\nStreaming trajectory for H-bonds, geometry, energy...")

hbond_counts = np.zeros(n_frames)
oh_sum = 0.0
oh_sum2 = 0.0
oh_count = 0
ang_sum = 0.0
ang_sum2 = 0.0
ang_count = 0
oh_hist_bins = np.linspace(0.8, 1.2, 200)
ang_hist_bins = np.linspace(90, 120, 200)
oh_hist = np.zeros(len(oh_hist_bins) - 1)
ang_hist = np.zeros(len(ang_hist_bins) - 1)

# Parse energies from comment lines (fast file scan)
energies = np.full(n_frames, np.nan)
n_atoms_file = len(u.atoms)
block = n_atoms_file + 2
try:
    with open(TRAJ, 'r') as f:
        raw = f.readlines()
    for fi in range(n_frames):
        line_idx = fi * STRIDE * block + 1
        if line_idx < len(raw):
            comment = raw[line_idx].strip()
            for token in comment.replace('=', ' ').split():
                try:
                    val = float(token)
                    if abs(val) > 1.0:
                        energies[fi] = val
                        break
                except ValueError:
                    pass
    valid_e = np.sum(~np.isnan(energies))
    print(f"  Energies found  : {valid_e} frames")
except Exception as ex:
    print(f"  Energy parse skipped: {ex}")

o_indices = O_sel.indices
h_indices = H_sel.indices

for fi, ts in enumerate(u.trajectory[::STRIDE]):
    pos = u.atoms.positions
    pos_O = pos[o_indices]
    pos_H = pos[h_indices]

    # ── O-H bond lengths + H-O-H angles (vectorized per molecule) ──
    O_pos = pos_O  # (N_W, 3)
    H1 = pos_H[0::2]  # (N_W, 3)
    H2 = pos_H[1::2]  # (N_W, 3)

    v1 = H1 - O_pos;
    v1 -= np.round(v1 / box3) * box3
    v2 = H2 - O_pos;
    v2 -= np.round(v2 / box3) * box3

    d1 = np.linalg.norm(v1, axis=1)  # (N_W,)
    d2 = np.linalg.norm(v2, axis=1)

    oh_vals = np.concatenate([d1, d2])
    oh_hist += np.histogram(oh_vals, bins=oh_hist_bins)[0]
    oh_sum += oh_vals.sum()
    oh_sum2 += (oh_vals ** 2).sum()
    oh_count += len(oh_vals)

    cos_a = np.sum(v1 * v2, axis=1) / (d1 * d2 + 1e-10)
    angles = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
    ang_hist += np.histogram(angles, bins=ang_hist_bins)[0]
    ang_sum += angles.sum()
    ang_sum2 += (angles ** 2).sum()
    ang_count += len(angles)

    # ── H-bonds (vectorized O-O) ──
    diff = pos_O[:, None, :] - pos_O[None, :, :]
    diff -= np.round(diff / box3) * box3
    oo_dist = np.sqrt((diff ** 2).sum(axis=-1))

    n_hb = 0
    for i in range(N_WATERS):
        cands = np.where((oo_dist[i] < HBOND_DIST) & (oo_dist[i] > 0.01))[0]
        if len(cands) == 0:
            continue
        for hj in [i * 2, i * 2 + 1]:
            oh_v = pos_H[hj] - pos_O[i]
            oh_v -= np.round(oh_v / box3) * box3
            oh_n = np.linalg.norm(oh_v)
            oo_vs = diff[i, cands]
            oo_ns = oo_dist[i, cands]
            cos_a = np.dot(oo_vs, oh_v) / (oo_ns * oh_n + 1e-10)
            n_hb += np.sum(np.degrees(np.arccos(np.clip(cos_a, -1, 1))) < HBOND_ANG)
    hbond_counts[fi] = n_hb / N_WATERS

oh_mean = oh_sum / oh_count
oh_std = np.sqrt(oh_sum2 / oh_count - oh_mean ** 2)
ang_mean = ang_sum / ang_count
ang_std = np.sqrt(ang_sum2 / ang_count - ang_mean ** 2)

print(f"  Mean O-H length       : {oh_mean:.4f} ± {oh_std:.4f} Å  (TIP3P: 0.9572 Å)")
print(f"  Mean H-O-H angle      : {ang_mean:.2f} ± {ang_std:.2f}°  (TIP3P: 104.52°)")
print(f"  Mean H-bonds/molecule : {hbond_counts.mean():.2f}  (ref: 3.5)")

# ─────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────
print("\nGenerating plots...")
plt.style.use('seaborn-v0_8-whitegrid')
fig = plt.figure(figsize=(16, 14))
fig.suptitle(f'CP2K PBE-D3 Water MD — 64 molecules, 12.4 Å box  (stride={STRIDE})',
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(3, 3, hspace=0.42, wspace=0.38)
c1, c2, c3 = '#1f77b4', '#d62728', '#2ca02c'


def ax_(pos, xlabel, ylabel, title, xlim=None):
    a = fig.add_subplot(pos)
    a.set(xlabel=xlabel, ylabel=ylabel, title=title)
    if xlim: a.set_xlim(xlim)
    return a


# RDF O-O
ax = ax_(gs[0, 0], 'r (Å)', 'g(r)', 'O-O RDF', (0, RDF_RMAX))
ax.plot(r_OO, g_OO, color=c1, lw=2)
ax.axvline(2.8, color='gray', ls='--', lw=1, label='Exp 2.80 Å')
ax.axvline(r_OO[peak_idx], color=c1, ls=':', lw=1.5, label=f'Peak {r_OO[peak_idx]:.2f} Å')
ax.legend(fontsize=8)

# RDF O-H
ax = ax_(gs[0, 1], 'r (Å)', 'g(r)', 'O-H RDF', (0, RDF_RMAX))
ax.plot(r_OH, g_OH, color=c2, lw=2)
ax.axvline(1.85, color='gray', ls='--', lw=1, label='Exp 1.85 Å')
ax.legend(fontsize=8)

# RDF H-H
ax = ax_(gs[0, 2], 'r (Å)', 'g(r)', 'H-H RDF', (0, RDF_RMAX))
ax.plot(r_HH, g_HH, color=c3, lw=2)
ax.axvline(2.4, color='gray', ls='--', lw=1, label='Exp 2.40 Å')
ax.legend(fontsize=8)

# MSD
ax = ax_(gs[1, 0], 'Time (ps)', 'MSD (Å²)', 'Mean Square Displacement')
ax.plot(times, msd, color=c1, lw=1.5, label='MSD')
ax.plot(times[s:e], slope * times[s:e] + intercept, 'r--', lw=2, label=f'D={D_SI:.1e} m²/s')
ax.legend(fontsize=8)

# H-bonds vs time
ax = ax_(gs[1, 1], 'Time (ps)', 'H-bonds/molecule', 'Hydrogen Bonds')
ax.plot(times, hbond_counts, color=c2, lw=1, alpha=0.6)
ax.axhline(hbond_counts.mean(), color='k', ls='--', lw=2, label=f'Mean={hbond_counts.mean():.2f}')
ax.axhline(3.5, color='gray', ls=':', lw=1.5, label='Ref=3.5')
ax.legend(fontsize=8)

# Energy
ax = ax_(gs[1, 2], 'Time (ps)', 'Energy (Ha)', 'Potential Energy')
valid = ~np.isnan(energies)
if valid.sum() > 0:
    ax.plot(times[valid], energies[valid], color=c3, lw=1, alpha=0.8)
    ax.axhline(np.nanmean(energies), color='k', ls='--', lw=2,
               label=f'Mean={np.nanmean(energies):.4f} Ha')
    ax.legend(fontsize=8)
else:
    ax.text(0.5, 0.5, 'No energy data\nin XYZ comments',
            ha='center', va='center', transform=ax.transAxes, fontsize=10)

# O-H bond length histogram
ax = ax_(gs[2, 0], 'O-H length (Å)', 'Counts', 'O-H Bond Length')
oh_centers = 0.5 * (oh_hist_bins[:-1] + oh_hist_bins[1:])
ax.bar(oh_centers, oh_hist, width=oh_centers[1] - oh_centers[0], color=c1, alpha=0.8, edgecolor='none')
ax.axvline(oh_mean, color='k', ls='--', lw=2, label=f'Mean={oh_mean:.4f} Å')
ax.axvline(0.9572, color='gray', ls=':', lw=1.5, label='TIP3P 0.9572 Å')
ax.legend(fontsize=8)

# H-O-H angle histogram
ax = ax_(gs[2, 1], 'H-O-H angle (°)', 'Counts', 'H-O-H Angle')
ang_centers = 0.5 * (ang_hist_bins[:-1] + ang_hist_bins[1:])
ax.bar(ang_centers, ang_hist, width=ang_centers[1] - ang_centers[0], color=c2, alpha=0.8, edgecolor='none')
ax.axvline(ang_mean, color='k', ls='--', lw=2, label=f'Mean={ang_mean:.2f}°')
ax.axvline(104.52, color='gray', ls=':', lw=1.5, label='TIP3P 104.52°')
ax.legend(fontsize=8)

# H-bond histogram
ax = ax_(gs[2, 2], 'H-bonds/molecule', 'Density', 'H-bond Distribution')
ax.hist(hbond_counts, bins=30, color=c3, alpha=0.8, density=True, edgecolor='none')
ax.axvline(hbond_counts.mean(), color='k', ls='--', lw=2, label=f'Mean={hbond_counts.mean():.2f}')
ax.legend(fontsize=8)

plt.savefig('water_analysis.png', dpi=150, bbox_inches='tight')
plt.show()
print("Plot saved: water_analysis.png")

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Frames analysed       : {n_frames}  (stride={STRIDE})")
print(f"Total time            : {times[-1]:.2f} ps")
print(f"O-O first peak        : {r_OO[peak_idx]:.3f} Å  (exp: 2.80 Å)")
print(f"O-O first min         : {r_OO[min_idx]:.3f} Å  (exp: 3.50 Å)")
print(f"Coordination number   : {coord_num:.2f}       (exp: 4.5)")
print(f"Diffusion coeff       : {D_SI:.2e} m²/s  (exp: 2.3e-9)")
print(f"Mean H-bonds/molecule : {hbond_counts.mean():.2f}       (exp: 3.5)")
print(f"Mean O-H length       : {oh_mean:.4f} ± {oh_std:.4f} Å")
print(f"Mean H-O-H angle      : {ang_mean:.2f} ± {ang_std:.2f}°")
if valid.sum() > 0:
    print(f"Mean energy           : {np.nanmean(energies):.6f} Ha")
print("=" * 60)
