from ase.calculators.calculator import Calculator, PropertyNotImplementedError, all_changes

from deepmd.infer import DeepPot

from pathlib import Path

import numpy as np

from typing import TYPE_CHECKING, ClassVar, Optional, Union, Callable

if TYPE_CHECKING:
    from ase import Atoms

__all__ = ["DP"]

class DP(Calculator):
    """Implementation of ASE deepmd calculator for polaron dynamics.

    Implemented properties are `energy`, `free_energy`, `forces`, `stress`

    Parameters

    ----------
    model_ener_force : Union[str, Path]
        path to the energy and forces model

    model_spin : Union[str, Path]
        path to the spin model

    label : str, optional

        calculator label, by default "DP"
    type_dict : dict[str, int], optional

        mapping of element types and their numbers, best left None and the calculator

        will infer this information from model, by default None

    neighbor_list : ase.neighborlist.NeighborList, optional

        The neighbor list object. If None, then build the native neighbor list.
    initial_charge_state : np.ndarray, optional

        Initial charge state of the system, by default None

    identifier : Callable, optional

        Callable to identify the charge state, by default None

    """

    name = "DP"
    implemented_properties: ClassVar[list[str]] = [
        "energy",
        "free_energy",
        "forces",
        "virial",
        "stress"
    ]

    def __init__(
            self,
            model_ener_force: Union[str, Path],
            model_spin: Union[str, Path],
            label: str = "DP",
            type_dict: Optional[dict[str, int]] = None,
            neighbor_list=None,
            initial_charge_state: Optional[np.ndarray] = None,
            identifier: Optional[Callable] = None,
            **kwargs,
    ) -> None:

        # 1. initialize the ener_force model and spin model

        Calculator.__init__(self, label=label, **kwargs)
        self.dp_ener_force = DeepPot(
            str(Path(model_ener_force).resolve()),
            neighbor_list=neighbor_list,
        )
        self.dp_spin = DeepPot(
            str(Path(model_spin).resolve()),
            neighbor_list=neighbor_list,
        )

        # 2. initialize the type_map

        if type_dict:
            self.type_dict = type_dict

        else:
            self.type_dict = dict(
                zip(self.dp_ener_force.get_type_map(), range(self.dp_ener_force.get_ntypes()))
            )

        # 3. set charge_state & spin

        self.charge_state = initial_charge_state  # record current charge state

        self.charge_state_history = [initial_charge_state]
        self.spin = None  # record current spin moment

        self.spin_history = [None]
        self.identifier = identifier

    def calculate(
            self,
            atoms: Optional["Atoms"] = None,
            properties: list[str] = ["energy", "forces", "virial"],
            system_changes: list[str] = all_changes,
    ) -> None:
        """Run calculation with deepmd models.

        Parameters

        ----------
        atoms : Optional[Atoms], optional

            atoms object to run the calculation on, by default None

        properties : list[str], optional

            unused, only for function signature compatibility,
            by default ["energy", "forces", "stress"]
        system_changes : list[str], optional

            unused, only for function signature compatibility, by default all_changes

        """

        # 0. prepare: coord, cell, symbols, atype

        if atoms is not None:
            self.atoms = atoms.copy()

        coord = self.atoms.get_positions().reshape([1, -1])
        if sum(self.atoms.get_pbc()) > 0:
            cell = self.atoms.get_cell().reshape([1, -1])
        else:
            cell = None

        symbols = self.atoms.get_chemical_symbols()
        atype = [self.type_dict[k] for k in symbols]

        # 1. get energy, forces, virial from the ener_force model

        e, f, v, _, _ = self.dp_ener_force.eval(coords=coord, cells=cell, atom_types=atype, atomic=True,
                                                aparam=self.charge_state)
        self.results["energy"] = e[0][0]
        self.results["free_energy"] = e[0][0]
        self.results["forces"] = f[0]
        self.results["virial"] = v[0].reshape(3, 3)

        # 2. get updated spin from the spin model

        _, _, _, spin, _ = self.dp_spin.eval(coords=coord, cells=cell, atom_types=atype, atomic=True,
                                             aparam=self.charge_state)
        self.spin = spin.flatten()
        self.spin_history.append(self.spin)
        self.charge_state = self.identifier(self.atoms, self.spin)
        self.charge_state_history.append(self.charge_state)

        # convert virial into stress for lattice relaxation

        if "stress" in properties:
            if sum(atoms.get_pbc()) > 0:
                # the usual convention (tensile stress is positive)
                # stress = -virial / volume

                stress = -0.5 * (v[0].copy() + v[0].copy().T) / atoms.get_volume()
                # Voigt notation

                self.results["stress"] = stress.flat[[0, 4, 8, 5, 2, 1]]
            else:
                raise PropertyNotImplementedError