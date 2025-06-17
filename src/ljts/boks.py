import random
import math
from src.ljts.Molecule import Molecule
import os

class Box:
    def __init__(self, length_x, length_y, length_z, density_liquid, density_vapor):
        self._lx = length_x  
        self._ly = length_y  
        self._lz = length_z  
        self._density_liquid = density_liquid  
        self._density_vapor = density_vapor 
        self._molecules = []  # list of molecules
        self._fill_box()  # populate the box with molecules
        self._last_energy = None  # cached potential energy

    def _fill_box(self):
        # Define region limits in the y-direction
        y_bottom = 0
        y_liquid_start = self._ly * 0.4
        y_liquid_end = self._ly * 0.6
        y_top = self._ly

        # Calculate volumes for liquid and vapor regions
        volume_liquid = self._lx * (y_liquid_end - y_liquid_start) * self._lz
        volume_vapor_single = self._lx * (y_liquid_start - y_bottom) * self._lz

        # Compute number of molecules to place in each region
        num_liquid = round(self._density_liquid * volume_liquid)
        num_vapor_total = round(self._density_vapor * volume_vapor_single * 2)

        # Add molecules to each region
        self._molecules += self._generate_molecules(num_vapor_total // 2, y_bottom, y_liquid_start)
        self._molecules += self._generate_molecules(num_liquid, y_liquid_start, y_liquid_end)
        self._molecules += self._generate_molecules(num_vapor_total // 2, y_liquid_end, y_top)

    def _generate_molecules(self, count, y_min, y_max):
        # Generate molecules randomly within a region
        new_molecules = []
        for _ in range(count):
            x = random.uniform(0, self._lx)
            y = random.uniform(y_min, y_max)
            z = random.uniform(0, self._lz)
            new_molecules.append(Molecule(x, y, z))
        return new_molecules

    def get_molecule_count(self):
        return len(self._molecules)

    def compute_potential_energy(self):
        """
        Calculate the total potential energy of the system
        """
        total_energy = 0.0
        N = len(self._molecules)
        for i in range(N):
            total_energy += 0.5 * self._energy_of_molecule(i)  # avoid double counting
        self._last_energy = total_energy
        return total_energy

    def get_last_energy(self):
        return self._last_energy

    def _energy_of_molecule(self, idx):
        cutoff = 2.5
        cutoff_squared = cutoff ** 2
        shift_value = 0.01631689

        mol_i = self._molecules[idx]
        energy = 0.0
        for j, mol_j in enumerate(self._molecules):
            if j == idx:
                continue
            dx = mol_i._x - mol_j._x
            dy = mol_i._y - mol_j._y
            dz = mol_i._z - mol_j._z

            dx -= round(dx / self._lx) * self._lx
            dy -= round(dy / self._ly) * self._ly
            dz -= round(dz / self._lz) * self._lz

            r_squared = dx*dx + dy*dy + dz*dz
            if r_squared < cutoff_squared and r_squared > 1e-12:
                inv_r2 = 1.0 / r_squared
                inv_r6 = inv_r2 ** 3
                inv_r12 = inv_r6 ** 2
                energy += 4 * (inv_r12 - inv_r6) + shift_value
        return energy

    def compute_distorted_energy(self, sx, sy, sz):

        total_energy = 0.0
        cutoff = 2.5
        cutoff_squared = cutoff ** 2
        shift_value = 0.01631689
        N = len(self._molecules)

        for i in range(N):
            for j in range(i+1, N):
                mol_i = self._molecules[i]
                mol_j = self._molecules[j]

                dx = (mol_i._x - mol_j._x) * sx
                dy = (mol_i._y - mol_j._y) * sy
                dz = (mol_i._z - mol_j._z) * sz

                #conditions on scaled box
                dx -= round(dx / (self._lx * sx)) * (self._lx * sx)
                dy -= round(dy / (self._ly * sy)) * (self._ly * sy)
                dz -= round(dz / (self._lz * sz)) * (self._lz * sz)

                r_squared = dx*dx + dy*dy + dz*dz
                if r_squared < cutoff_squared and r_squared > 1e-12:
                    inv_r2 = 1.0 / r_squared
                    inv_r6 = inv_r2 ** 3
                    inv_r12 = inv_r6 ** 2
                    total_energy += 4 * (inv_r12 - inv_r6) + shift_value
        return total_energy

    def get_surface_area(self):
        
        return 2 * self._lx * self._lz

    def mc_sweep(self, b, T):
        # Perform one MC sweep (N attempted moves)
        accepted = 0
        for _ in range(len(self._molecules)):
            if self.mc_move(b, T):
                accepted += 1
        return accepted

    def mc_move(self, b, T):
        # single Monte Carlo move
        i = random.randrange(len(self._molecules))
        mol = self._molecules[i]
        old_x, old_y, old_z = mol._x, mol._y, mol._z

        old_energy = self._energy_of_molecule(i)

        # Random displacement
        mol._x = (mol._x + random.uniform(-b, b)) % self._lx
        mol._y = (mol._y + random.uniform(-b, b)) % self._ly
        mol._z = (mol._z + random.uniform(-b, b)) % self._lz

        new_energy = self._energy_of_molecule(i)
        delta_E = new_energy - old_energy

        if delta_E <= 0 or random.random() < math.exp(-delta_E / T):
            return True  # Accept move
        else:
            mol._x, mol._y, mol._z = old_x, old_y, old_z  # Revert move
            return False

    def export_xyz(self, filepath, comment=""):
        # Export current configuration to an XYZ file
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"{len(self._molecules)}\n")
            f.write(comment + "\n")
            for mol in self._molecules:
                x, y, z = mol.get_position()
                f.write(f"Position x:{x:.6f} y:{y:.6f} z:{z:.6f}\n")

    def append_xyz_frame(self, filepath, comment=""):
        # Append current frame to trajectory in XYZ format
        with open(filepath, 'a', encoding='utf-8') as f:
            f.write(f"{len(self._molecules)}\n")
            f.write(comment + "\n")
            for mol in self._molecules:
                x, y, z = mol.get_position()
                f.write(f"H {x:.6f} {y:.6f} {z:.6f}\n") #for ovito visiulization