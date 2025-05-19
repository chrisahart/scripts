# SPDX-License-Identifier: LGPL-3.0-or-later
# A modified version of original DP ASE calculator for polaron dynamics
"""ASE calculator interface module."""
import numpy as np
from pathlib import (
    Path,
)
from typing import (
    TYPE_CHECKING,
    ClassVar,
    Optional,
    Union,
    Callable
)

from ase.calculators.calculator import (
    Calculator,
    PropertyNotImplementedError,
    all_changes,
)

from deepmd.infer import (
    DeepPot,
)

if TYPE_CHECKING:
    from ase import (
        Atoms,
    )

__all__ = ["DP"]


class DP(Calculator):
    """Implementation of ASE deepmd calculator.

    Implemented properties are `energy`, `forces` and `stress`

    Parameters
    ----------
    model : Union[str, Path]
        path to the model
    label : str, optional
        calculator label, by default "DP"
    type_dict : dict[str, int], optional
        mapping of element types and their numbers, best left None and the calculator
        will infer this information from model, by default None
    neighbor_list : ase.neighborlist.NeighborList, optional
        The neighbor list object. If None, then build the native neighbor list.
    head : Union[str, None], optional
        a specific model branch choosing from pretrained model, by default None

    Examples
    --------
    Compute potential energy

    >>> from ase import Atoms
    >>> from deepmd.tf.calculator import DP
    >>> water = Atoms('H2O',
    >>>             positions=[(0.7601, 1.9270, 1),
    >>>                        (1.9575, 1, 1),
    >>>                        (1., 1., 1.)],
    >>>             cell=[100, 100, 100],
    >>>             calculator=DP(model="frozen_model.pb"))
    >>> print(water.get_potential_energy())
    >>> print(water.get_forces())

    Run BFGS structure optimization

    >>> from ase.optimize import BFGS
    >>> dyn = BFGS(water)
    >>> dyn.run(fmax=1e-6)
    >>> print(water.get_positions())
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
            model: Union[str, "Path"],
            label: str = "DP",
            type_dict: Optional[dict[str, int]] = None,
            neighbor_list=None,
            force_head: str = None,
            charge_head: str = None,
            initial_charge_state: np.ndarray = None,
            identifier: Optional[Callable] = None,
            **kwargs,
    ) -> None:

        # 1. initialize the ener_force model, and charge model
        Calculator.__init__(self, label=label, **kwargs)
        self.dp_force = DeepPot(
            str(Path(model).resolve()),
            neighbor_list=neighbor_list,
            head=force_head,
        )
        self.dp_charge = DeepPot(
            str(Path(model).resolve()),
            neighbor_list=neighbor_list,
            head=charge_head,
        )

        # 2. initialize the type_map
        if type_dict:
            self.type_dict = type_dict
        else:
            self.type_dict = dict(
                zip(self.dp_force.get_type_map(), range(self.dp_force.get_ntypes()))
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
        """Run calculation with deepmd model.

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

        # 1. get energy, forces, virial
        e, f, v, _, _ = self.dp_force.eval(coords=coord, cells=cell, atom_types=atype, atomic=True,
                                           aparam=self.charge_state)
        self.results["energy"] = e[0][0]
        # see https://gitlab.com/ase/ase/-/merge_requests/2485
        self.results["free_energy"] = e[0][0]
        self.results["forces"] = f[0]
        self.results["virial"] = v[0].reshape(3, 3)

        # 2. get updated charge state
        _, _, _, spin, _ = self.dp_charge.eval(coords=coord, cells=cell, atom_types=atype, atomic=True,
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