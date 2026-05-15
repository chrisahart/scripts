import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.transformations import set_dimensions
from MDAnalysis.lib.distances import apply_PBC
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

# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/cubic/64/md/pbe/md/from-tip3p-1-ns'
# TRAJ = '{}/water-pos-1.xyz'.format(folder)
# box_cubic = 12.4
# BOX = [box_cubic, box_cubic, box_cubic, 90.0, 90.0, 90.0]
# TRAJ_START = 50000
# RDF_RMAX = box_cubic / 2
# N_WATERS = 64

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/cubic/128/md/pbe/from-tip3p-1-ns'
TRAJ = '{}/water-pos-1.xyz'.format(folder)
box_cubic = 15.65
BOX = [box_cubic, box_cubic, box_cubic, 90.0, 90.0, 90.0]
TRAJ_START = 50000
RDF_RMAX = box_cubic / 2
N_WATERS = 128

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/pbe-d3-bj-c9-pbe/from-tip3p-1-ns'
TRAJ = '{}/water-pos-1.xyz'.format(folder)
box_cubic = 15.65
BOX = [box_cubic, box_cubic, box_cubic, 90.0, 90.0, 90.0]
TRAJ_START = 50000
RDF_RMAX = box_cubic / 2
N_WATERS = 128

STRIDE = 1
RDF_BINS = 200
save_fig = True
save_fig = False
OUTDIR = folder
box3 = np.array(BOX[:3])

# Water molecular mass: O + 2H (g/mol)
MW_WATER = 15.9994 + 2 * 1.0079  # 18.0152 g/mol
AVOGADRO = 6.02214076e23

print("=" * 60)
print("CP2K Water Box MD Analysis")
print(f"Output folder  : {OUTDIR}")
print(f"TRAJ_START     : frame {TRAJ_START}  (equilibration skipped)")
print("=" * 60)

# ─────────────────────────────────────────────
# LOAD
# ─────────────────────────────────────────────
u = mda.Universe(TRAJ, format='XYZ')
u.trajectory.add_transformations(set_dimensions(BOX))

n_frames_total = len(u.trajectory)
n_frames = (n_frames_total - TRAJ_START) // STRIDE

print(f"Total frames   : {n_frames_total}")
print(f"Skipped        : {TRAJ_START} frames (equilibration)")
print(f"Stride         : {STRIDE}")
print(f"Frames used    : {n_frames}")

O_sel = u.select_atoms('name O')
H_sel = u.select_atoms('name H')
prod_slice = slice(TRAJ_START, n_frames_total, STRIDE)

# ─────────────────────────────────────────────
# 1. DENSITY (NVT — computed once from fixed box)
# ─────────────────────────────────────────────
print("\nCalculating density (NVT — fixed box)...")

mass_g = (N_WATERS * MW_WATER) / AVOGADRO  # grams
vol_A3 = BOX[0] * BOX[1] * BOX[2]  # Å³
vol_cm3 = vol_A3 * 1e-24  # cm³
rho = mass_g / vol_cm3  # g/cm³
rho_ref = 1.0
rho_err = abs(rho - rho_ref) / rho_ref * 100

print(f"  Box volume  : {vol_A3:.4f} Å³")
print(f"  Mass        : {mass_g * AVOGADRO / N_WATERS:.4f} g/mol per molecule")
print(f"  Density     : {rho:.4f} g/cm³")
print(f"  Reference   : {rho_ref:.4f} g/cm³  (liquid water 298 K)")
print(f"  Deviation   : {rho_err:.2f}%  {'✅' if rho_err < 5 else '⚠️'}")

# ─────────────────────────────────────────────
# 2. RDFs
# ─────────────────────────────────────────────
print("\nCalculating RDFs...")

rdf_OO = InterRDF(O_sel, O_sel, nbins=RDF_BINS, range=(0.5, RDF_RMAX), exclusion_block=(1, 1))
rdf_OH = InterRDF(O_sel, H_sel, nbins=RDF_BINS, range=(0.5, RDF_RMAX))
rdf_HH = InterRDF(H_sel, H_sel, nbins=RDF_BINS, range=(0.5, RDF_RMAX), exclusion_block=(1, 1))

rdf_OO.run(start=TRAJ_START, step=STRIDE, verbose=True)
rdf_OH.run(start=TRAJ_START, step=STRIDE, verbose=True)
rdf_HH.run(start=TRAJ_START, step=STRIDE, verbose=True)

r_OO, g_OO = rdf_OO.results.bins, rdf_OO.results.rdf
r_OH, g_OH = rdf_OH.results.bins, rdf_OH.results.rdf
r_HH, g_HH = rdf_HH.results.bins, rdf_HH.results.rdf

dr = r_OO[1] - r_OO[0]
peak_OO = r_OO[int(2.0 / dr) + np.argmax(g_OO[int(2.0 / dr):int(4.0 / dr)])]
peak_OH = r_OH[int(0.8 / dr) + np.argmax(g_OH[int(0.8 / dr):int(2.5 / dr)])]
peak_HH = r_HH[int(1.5 / dr) + np.argmax(g_HH[int(1.5 / dr):int(3.0 / dr)])]

# Coordination number — integrate g_OO to first minimum
peak_OO_idx = int(2.0 / dr) + np.argmax(g_OO[int(2.0 / dr):int(4.0 / dr)])
min_OO_idx = peak_OO_idx + np.argmin(g_OO[peak_OO_idx:peak_OO_idx + int(2.0 / dr)])
rho_O = len(O_sel) / box3.prod()
try:
    coord_num = 4 * np.pi * rho_O * np.trapezoid(g_OO[:min_OO_idx] * r_OO[:min_OO_idx] ** 2, r_OO[:min_OO_idx])
except AttributeError:
    coord_num = 4 * np.pi * rho_O * np.trapz(g_OO[:min_OO_idx] * r_OO[:min_OO_idx] ** 2, r_OO[:min_OO_idx])

print(f"  O-O first peak : {peak_OO:.3f} Å  (ref: 2.80 Å)")
print(f"  O-H first peak : {peak_OH:.3f} Å  (ref: 1.85 Å)")
print(f"  H-H first peak : {peak_HH:.3f} Å  (ref: 2.40 Å)")
print(f"  Coordination N : {coord_num:.2f}  (ref: 4.5)")

# ─────────────────────────────────────────────
# 3. O-H BOND LENGTHS + H-O-H ANGLES
#
# PBC strategy:
#   Step 1 — apply_PBC(pos, dims): wraps absolute coordinates into
#            the primary box [0, box]. Ensures O and H positions are
#            in the same image before computing displacements.
#   Step 2 — MIC on displacement vectors: v -= round(v/box)*box
#            maps each O→H vector into [-box/2, +box/2], giving the
#            shortest image distance. This is the correct minimum
#            image convention for bond vectors.
#   Note: for intramolecular bonds (~0.97 Å << box/2 ~6.2 Å),
#         step 1 is redundant but harmless. Step 2 alone is sufficient.
# ─────────────────────────────────────────────
print("\nCalculating O-H lengths and H-O-H angles...")

oh_hist_bins = np.linspace(0.80, 1.20, 300)
ang_hist_bins = np.linspace(90.0, 125.0, 300)
oh_hist = np.zeros(len(oh_hist_bins) - 1)
ang_hist = np.zeros(len(ang_hist_bins) - 1)
oh_sum = 0.0;
oh_sum2 = 0.0;
oh_n = 0
ang_sum = 0.0;
ang_sum2 = 0.0;
ang_n = 0

o_idx = O_sel.indices
h_idx = H_sel.indices

for ts in u.trajectory[prod_slice]:
    dims = ts.dimensions
    pos = u.atoms.positions.copy()

    # Step 1: wrap coordinates into primary box
    pos = apply_PBC(pos, dims)
    pos_O = pos[o_idx]
    pos_H = pos[h_idx]

    H1 = pos_H[0::2]
    H2 = pos_H[1::2]

    # Step 2: MIC on displacement vectors
    v1 = H1 - pos_O;
    v1 -= np.round(v1 / box3) * box3
    v2 = H2 - pos_O;
    v2 -= np.round(v2 / box3) * box3

    d1 = np.linalg.norm(v1, axis=1)
    d2 = np.linalg.norm(v2, axis=1)
    oh_vals = np.concatenate([d1, d2])

    oh_hist += np.histogram(oh_vals, bins=oh_hist_bins)[0]
    oh_sum += oh_vals.sum();
    oh_sum2 += (oh_vals ** 2).sum();
    oh_n += len(oh_vals)

    cos_a = np.sum(v1 * v2, axis=1) / (d1 * d2 + 1e-10)
    angles = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))

    ang_hist += np.histogram(angles, bins=ang_hist_bins)[0]
    ang_sum += angles.sum();
    ang_sum2 += (angles ** 2).sum();
    ang_n += len(angles)

oh_mean = oh_sum / oh_n
oh_std = np.sqrt(oh_sum2 / oh_n - oh_mean ** 2)
ang_mean = ang_sum / ang_n
ang_std = np.sqrt(ang_sum2 / ang_n - ang_mean ** 2)

print(f"  Mean O-H length  : {oh_mean:.4f} ± {oh_std:.4f} Å")
print(f"  Mean H-O-H angle : {ang_mean:.2f} ± {ang_std:.2f}°")

# ─────────────────────────────────────────────
# FIGURE 1: O-O RDF
# ─────────────────────────────────────────────
fig_rdf_oo, ax_rdf_oo = plt.subplots(figsize=(7, 5))
ax_rdf_oo.plot(r_OO, g_OO, 'k-')
ax_rdf_oo.axhline(1.0, color='lightgray', ls='-', lw=1, zorder=0)
ax_rdf_oo.set_xlabel(r'r$_{\mathrm{OO}}$ (Å)')
ax_rdf_oo.set_ylabel(r'g$_{\mathrm{OO}}$(r)')
ax_rdf_oo.set_xlim([0, RDF_RMAX])
ax_rdf_oo.set_ylim([0, 4])
fig_rdf_oo.tight_layout()
if save_fig: fig_rdf_oo.savefig('{}/rdf_oo.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# FIGURE 2: O-H RDF
# ─────────────────────────────────────────────
fig_rdf_oh, ax_rdf_oh = plt.subplots(figsize=(7, 5))
ax_rdf_oh.plot(r_OH, g_OH, 'k-')
ax_rdf_oh.axhline(1.0, color='lightgray', ls='-', lw=1, zorder=0)
ax_rdf_oh.set_xlabel(r'r$_{\mathrm{OH}}$ (Å)')
ax_rdf_oh.set_ylabel(r'g$_{\mathrm{OH}}$(r)')
ax_rdf_oh.set_xlim([0, RDF_RMAX])
ax_rdf_oh.set_ylim([0, 2.15])
fig_rdf_oh.tight_layout()
if save_fig: fig_rdf_oh.savefig('{}/rdf_oh.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# FIGURE 3: H-H RDF
# ─────────────────────────────────────────────
fig_rdf_hh, ax_rdf_hh = plt.subplots(figsize=(7, 5))
ax_rdf_hh.plot(r_HH, g_HH, 'k-')
ax_rdf_hh.axhline(1.0, color='lightgray', ls='-', lw=1, zorder=0)
ax_rdf_hh.set_xlabel(r'r$_{\mathrm{HH}}$ (Å)')
ax_rdf_hh.set_ylabel(r'g$_{\mathrm{HH}}$(r)')
ax_rdf_hh.set_xlim([0, RDF_RMAX])
ax_rdf_hh.set_ylim([0, 3.15])
fig_rdf_hh.tight_layout()
if save_fig: fig_rdf_hh.savefig('{}/rdf_hh.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# FIGURE 4: O-H BOND LENGTH
# ─────────────────────────────────────────────
oh_centers = 0.5 * (oh_hist_bins[:-1] + oh_hist_bins[1:])
fig_oh_bond, ax_oh_bond = plt.subplots(figsize=(7, 5))
ax_oh_bond.bar(oh_centers, oh_hist, width=oh_centers[1] - oh_centers[0],
               color='black', edgecolor='none')
ax_oh_bond.axvline(oh_mean, color='gray', ls='-', lw=1, alpha=0.85,
                   label=f'Mean = {oh_mean:.4f} Å')
ax_oh_bond.set_xlabel('O-H Bond Length (Å)')
ax_oh_bond.set_ylabel('Counts')
ax_oh_bond.set_xlim([0.8, 1.20])
fig_oh_bond.tight_layout()
if save_fig: fig_oh_bond.savefig('{}/oh_bond.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# FIGURE 5: H-O-H ANGLE
# ─────────────────────────────────────────────
ang_centers = 0.5 * (ang_hist_bins[:-1] + ang_hist_bins[1:])
fig_hoh_angle, ax_hoh_angle = plt.subplots(figsize=(7, 5))
ax_hoh_angle.bar(ang_centers, ang_hist, width=ang_centers[1] - ang_centers[0],
                 color='black', edgecolor='none')
ax_hoh_angle.axvline(ang_mean, color='gray', ls='-', lw=1, alpha=0.5,
                     label=f'Mean = {ang_mean:.2f}°')
ax_hoh_angle.set_xlabel('H-O-H Angle (°)')
ax_hoh_angle.set_ylabel('Counts')
ax_hoh_angle.set_xlim([90, 125])
fig_hoh_angle.tight_layout()
if save_fig: fig_hoh_angle.savefig('{}/hoh_angle.png'.format(folder), dpi=300)

# ─────────────────────────────────────────────
# FIGURE 6: COMBINED RDF — 1x3 stacked, shared x-axis
# ─────────────────────────────────────────────
fig_rdf_all, (ax_oo, ax_oh, ax_hh) = plt.subplots(
    3, 1, figsize=(7, 10),
    sharex=True,
    gridspec_kw={"hspace": 0.00}
)

ax_oo.plot(r_OO, g_OO, "k-")
ax_oo.axhline(1.0, color="lightgray", ls="-", lw=1, zorder=0)
ax_oo.set_ylabel(r"g$_{\mathrm{OO}}$(r)")
ax_oo.set_xlim([1, 8])
ax_oo.set_ylim([0, 4])

ax_oh.plot(r_OH, g_OH, "k-")
ax_oh.axhline(1.0, color="lightgray", ls="-", lw=1, zorder=0)
ax_oh.set_ylabel(r"g$_{\mathrm{OH}}$(r)")
ax_oh.set_xlim([1, 8])
ax_oh.set_ylim([0, 3.99])

ax_hh.plot(r_HH, g_HH, "k-")
ax_hh.axhline(1.0, color="lightgray", ls="-", lw=1, zorder=0)
ax_hh.set_ylabel(r"g$_{\mathrm{HH}}$(r)")
ax_hh.set_xlim([1, 8])
ax_hh.set_ylim([0, 3.99])
ax_hh.set_xlabel(r"r (Å)")

fig_rdf_all.tight_layout()
if save_fig: fig_rdf_all.savefig("{}/rdf_all.png".format(folder), dpi=300)

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Total frames     : {n_frames_total}")
print(f"Equilibration    : frames 0 → {TRAJ_START - 1} (skipped)")
print(f"Production       : frames {TRAJ_START} → {n_frames_total - 1}  ({n_frames} used)")
print(f"Density          : {rho:.4f} g/cm³  (ref: 1.0000, Δ={rho_err:.2f}%)")
print(f"O-O first peak   : {peak_OO:.3f} Å  (exp: 2.80 Å)")
print(f"O-H first peak   : {peak_OH:.3f} Å  (exp: 1.85 Å)")
print(f"H-H first peak   : {peak_HH:.3f} Å  (exp: 2.40 Å)")
print(f"Coordination N   : {coord_num:.2f}       (exp: 4.5)")
print(f"Mean O-H length  : {oh_mean:.4f} ± {oh_std:.4f} Å")
print(f"Mean H-O-H angle : {ang_mean:.2f} ± {ang_std:.2f}°")
print("=" * 60)

if __name__ == "__main__":
    print('Finished.')
    plt.show()