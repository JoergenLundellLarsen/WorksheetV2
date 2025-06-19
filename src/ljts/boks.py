import random
import math
from src.ljts.Molecule import Molecule
import os

class Box:
    """
    Simulation box for a Lennard-Jones TS fluid system in slab geometry.

    Responsible for:
    - Initializing molecule positions in vapor-liquid-vapor regions
    - Performing Monte Carlo sweeps
    - Computing potential energy and surface tension via test-area method
    - Exporting simulation snapshots
    """

    def __init__(self, length_x, length_y, length_z, density_liquid, density_vapor):
        self._lx = length_x
        self._ly = length_y
        self._lz = length_z
        self._density_liquid = density_liquid
        self._density_vapor = density_vapor
        self._molecules = []
        self._fill_box()
        self._last_energy = None

    def _fill_box(self):
        """
        Fill the box with molecules divided into vapor-liquid-vapor slabs along y-axis.
        """
        y_bottom = 0
        y_liquid_start = self._ly * 0.4
        y_liquid_end = self._ly * 0.6
        y_top = self._ly

        # Compute volumes of regions
        volume_liquid = self._lx * (y_liquid_end - y_liquid_start) * self._lz
        volume_vapor_single = self._lx * (y_liquid_start - y_bottom) * self._lz

        # Number of molecules from density
        num_liquid = round(self._density_liquid * volume_liquid)
        num_vapor_total = round(self._density_vapor * volume_vapor_single * 2)

        # Distribute molecules to each region
        self._molecules += self._generate_molecules(num_vapor_total // 2, y_bottom, y_liquid_start)
        self._molecules += self._generate_molecules(num_liquid, y_liquid_start, y_liquid_end)
        self._molecules += self._generate_molecules(num_vapor_total // 2, y_liquid_end, y_top)

    def _generate_molecules(self, count, y_min, y_max):
        """
        Randomly generate molecules in a subregion.
        """
        return [
            Molecule(
                random.uniform(0, self._lx),
                random.uniform(y_min, y_max),
                random.uniform(0, self._lz)
            ) for _ in range(count)
        ]

    def get_molecule_count(self):
        """Return total number of molecules."""
        return len(self._molecules)

    def compute_potential_energy(self):
        """
        Compute total potential energy of the system using pairwise interactions.
        """
        total_energy = sum(0.5 * self._energy_of_molecule(i) for i in range(len(self._molecules)))
        self._last_energy = total_energy
        return total_energy

    def get_last_energy(self):
        """Return last computed potential energy (if cached)."""
        return self._last_energy

    def _energy_of_molecule(self, idx):
        """
        Compute pair interaction energy for one molecule with all others.
        """
        cutoff = 2.5
        cutoff_squared = cutoff ** 2
        shift = 0.01631689  # energy shift to make potential zero at cutoff
        mol_i = self._molecules[idx]
        energy = 0.0

        for j, mol_j in enumerate(self._molecules):
            if j == idx:
                continue

            # Distance with periodic boundary conditions
            dx = mol_i._x - mol_j._x
            dy = mol_i._y - mol_j._y
            dz = mol_i._z - mol_j._z

            dx -= round(dx / self._lx) * self._lx
            dy -= round(dy / self._ly) * self._ly
            dz -= round(dz / self._lz) * self._lz

            r2 = dx*dx + dy*dy + dz*dz
            if 1e-12 < r2 < cutoff_squared:
                inv_r2 = 1.0 / r2
                inv_r6 = inv_r2 ** 3
                inv_r12 = inv_r6 ** 2
                energy += 4 * (inv_r12 - inv_r6) + shift
        return energy

    def compute_distorted_energy(self, sx, sy, sz):
        """
        Compute total energy of the system under a box deformation.
        Used in the test-area method.
        """
        total_energy = 0.0
        cutoff_squared = 2.5 ** 2
        shift = 0.01631689

        for i in range(len(self._molecules)):
            for j in range(i + 1, len(self._molecules)):
                mol_i = self._molecules[i]
                mol_j = self._molecules[j]

                # Scale coordinates
                dx = (mol_i._x - mol_j._x) * sx
                dy = (mol_i._y - mol_j._y) * sy
                dz = (mol_i._z - mol_j._z) * sz

                # Apply periodic conditions in scaled box
                dx -= round(dx / (self._lx * sx)) * (self._lx * sx)
                dy -= round(dy / (self._ly * sy)) * (self._ly * sy)
                dz -= round(dz / (self._lz * sz)) * (self._lz * sz)

                r2 = dx*dx + dy*dy + dz*dz
                if 1e-12 < r2 < cutoff_squared:
                    inv_r2 = 1.0 / r2
                    inv_r6 = inv_r2 ** 3
                    inv_r12 = inv_r6 ** 2
                    total_energy += 4 * (inv_r12 - inv_r6) + shift
        return total_energy

    def get_surface_area(self):
        """
        Return total interfacial area (two surfaces).
        """
        return 2 * self._lx * self._lz

    def mc_sweep(self, b, T):
        """
        Perform one full Monte Carlo sweep.
        """
        accepted = 0
        for _ in range(len(self._molecules)):
            accepted += self.mc_move(b, T)
        return accepted

    def mc_move(self, b, T):
        """
        Attempt a random displacement of one molecule.
        Accept or reject using the Metropolis criterion.
        """
        i = random.randrange(len(self._molecules))
        mol = self._molecules[i]
        old_pos = mol._x, mol._y, mol._z
        old_E = self._energy_of_molecule(i)

        # Random trial move
        mol._x = (mol._x + random.uniform(-b, b)) % self._lx
        mol._y = (mol._y + random.uniform(-b, b)) % self._ly
        mol._z = (mol._z + random.uniform(-b, b)) % self._lz

        new_E = self._energy_of_molecule(i)
        dE = new_E - old_E

        # Accept or reject move
        if dE <= 0 or random.random() < math.exp(-dE / T):
            return 1
        else:
            mol._x, mol._y, mol._z = old_pos
            return 0

    def export_xyz(self, filepath, comment=""):
        """
        Write current configuration to an XYZ file.
        """
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"{len(self._molecules)}\n")
            f.write(comment + "\n")
            for mol in self._molecules:
                x, y, z = mol.get_position()
                f.write(f"Position x:{x:.6f} y:{y:.6f} z:{z:.6f}\n")

    def append_xyz_frame(self, filepath, comment=""):
        """
        Append a frame to the trajectory file (XYZ format).
        """
        with open(filepath, 'a', encoding='utf-8') as f:
            f.write(f"{len(self._molecules)}\n")
            f.write(comment + "\n")
            for mol in self._molecules:
                x, y, z = mol.get_position()
                f.write(f"H {x:.6f} {y:.6f} {z:.6f}\n")  # for OVITO visualization
