import ase
from ase.md.md import MolecularDynamics
from ase.md.nose_hoover_chain import NoseHooverChainNVT

nvt_dyn = NoseHooverChainNVT()


print(ase.__version__)