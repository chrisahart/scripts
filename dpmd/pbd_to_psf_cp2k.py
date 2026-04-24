folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/testing/xyz'
n_waters = 64
n_atoms = n_waters * 3
n_bonds = n_waters * 2

with open('{}/system.psf'.format(folder), 'w') as f:
    # Header
    f.write("PSF\n\n")
    f.write("       1 !NTITLE\n")
    f.write(" REMARKS TIP3P water box 64 molecules\n\n")

    # Atoms - strict column format for CP2K Fortran reader
    # I8, 1X, A4, 1X, I4, 1X, A4, 1X, A4, 1X, A4, 2X, F10.6, 2X, F8.4, 2X, I4
    f.write(f"{n_atoms:8d} !NATOM\n")
    atom_idx = 1
    for mol in range(1, n_waters + 1):
        f.write(f"{atom_idx:8d} {'WAT':<4s} {mol:<4d} {'WAT':<4s} {'O':<4s} {'O':<4s}  {-0.834000:10.6f}  {15.9994:8.4f}  {'0':>4}\n")
        f.write(f"{atom_idx+1:8d} {'WAT':<4s} {mol:<4d} {'WAT':<4s} {'H1':<4s} {'H':<4s}  { 0.417000:10.6f}  { 1.0079:8.4f}  {'0':>4}\n")
        f.write(f"{atom_idx+2:8d} {'WAT':<4s} {mol:<4d} {'WAT':<4s} {'H2':<4s} {'H':<4s}  { 0.417000:10.6f}  { 1.0079:8.4f}  {'0':>4}\n")
        atom_idx += 3

    # Bonds - 8 indices per line
    f.write(f"\n{n_bonds:8d} !NBOND: bonds\n")
    count = 0
    for mol in range(n_waters):
        o  = mol * 3 + 1
        h1 = mol * 3 + 2
        h2 = mol * 3 + 3
        f.write(f"{o:8d}{h1:8d}")
        count += 1
        if count % 4 == 0:
            f.write("\n")
        f.write(f"{o:8d}{h2:8d}")
        count += 1
        if count % 4 == 0:
            f.write("\n")

    f.write("\n\n       0 !NTHETA: angles\n")
    f.write("\n\n       0 !NPHI: dihedrals\n")
    f.write("\n\n       0 !NIMPHI: impropers\n")
    f.write("\n\n       0 !NDON: donors\n")
    f.write("\n\n       0 !NACC: acceptors\n")
    f.write("\n\n       0 !NNB\n")

print("PSF written successfully")