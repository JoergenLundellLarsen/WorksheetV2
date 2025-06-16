import random
import math
from src.ljts.Molecule import Molecule

class Box:

    def __init__(self, length_x, length_y, length_z, density_liquid, density_vapor):
    
        self._lx = length_x
        self._ly = length_y
        self._lz = length_z
        self._density_liquid = density_liquid
        self._density_vapor = density_vapor
        self._molecules = []

        # Fill the box with molecules in vapor and liquid regions
        self._fill_box()

    def _fill_box(self):
        """
        Fill the box with molecules, placing them randomly in three regions along the y axis:
        - Bottom 40%: vapor
        - Middle 20%: liquid
        - Top 40%: vapor
        doe
        """

        #calculate y boundaries for the three regions
        y_bottom = 0
        y_liquid_start = self._ly * 0.4
        y_liquid_end = self._ly * 0.6
        y_top = self._ly

        #calculate the volumes for the liquid and vapor regions
        volume_liquid = self._lx * (y_liquid_end - y_liquid_start) * self._lz
        volume_vapor_single = self._lx * (y_liquid_start - y_bottom) * self._lz  # Each vapor region

        # calculate number of molecules in each region based on destity
        num_liquid = round(self._density_liquid * volume_liquid)
        num_vapor_total = round(self._density_vapor * volume_vapor_single * 2)  # Total vapor molecules

        # Generate and add molecules for each region
        self._molecules += self._generate_molecules(num_vapor_total // 2, y_bottom, y_liquid_start)
        self._molecules += self._generate_molecules(num_liquid, y_liquid_start, y_liquid_end)
        self._molecules += self._generate_molecules(num_vapor_total // 2, y_liquid_end, y_top)

    def _generate_molecules(self, count, y_min, y_max):
        """
        Generate  molecules randomly positioned in the box between y_min and y_max.
        """
        new_molecules = []
        for _ in range(count):
            x = random.uniform(0, self._lx)
            y = random.uniform(y_min, y_max)
            z = random.uniform(0, self._lz)
            new_molecules.append(Molecule(x, y, z))
        return new_molecules

    def get_molecule_count(self):
        # Return the total number of molecules in the box
        return len(self._molecules)

    def compute_potential_energy(self):
        """
        Calculate the total potential energy of the system
        """
        total_energy = 0.0
        cutoff = 2.5  # Potential cutoff distance
        cutoff_squared = cutoff ** 2 
        shift_value = 0.01631689  # Energy shift so potential = 0 at cutoff

        num_molecules = len(self._molecules)

        # Double loop over all unique molecule pairs (i < j)
        for i in range(num_molecules):
            mol_i = self._molecules[i]

            for j in range(i + 1, num_molecules):
                mol_j = self._molecules[j]

                # Calculate distance between molecules, considering periodic boundary conditions
                dx = mol_i._x - mol_j._x
                dy = mol_i._y - mol_j._y
                dz = mol_i._z - mol_j._z

                # Minimum image convention (periodic boundary conditions)
                dx -= round(dx / self._lx) * self._lx
                dy -= round(dy / self._ly) * self._ly
                dz -= round(dz / self._lz) * self._lz

                # Compute squared distance
                r_squared = dx * dx + dy * dy + dz * dz

                if r_squared < cutoff_squared:
                    # Lennard-Jones potential: u_LJ(r) = 4 * (1/r^12 - 1/r^6) + shift_value
                    inv_r2 = 1.0 / r_squared
                    inv_r6 = inv_r2 ** 3
                    inv_r12 = inv_r6 ** 2

                    energy_ij = 4 * (inv_r12 - inv_r6) + shift_value
                    total_energy += energy_ij

        return total_energy

    def mc_sweep(self, b, T):
        """
        one Monte Carlo sweep:
        Try to move each molecule once
        Return the number of accepted moves
        """
        N = len(self._molecules)
        accepted = 0
        for _ in range(N):
            if self.mc_move(b, T):
                accepted += 1
        return accepted

    def mc_move(self, b, T):
        """
        Pick a random molecule
        Calculate its energy before the move (interactions with others)
        Move it randomly within [-b, b] in each dimension, using PBC
        Calculate its new energy
        Accept the move if energy decreased, or with exp(-ΔE/T) if energy increased
        If rejected, restore the old position
        Returns True if the move was accepted, False otherwise.
        """
        i = random.randint(0, len(self._molecules) - 1)
        mol = self._molecules[i]
        old_x, old_y, old_z = mol._x, mol._y, mol._z

        # Calculate energy before the move
        old_energy = self._energy_of_molecule(i)

        # Make a trial move (add random displacement in each direction)
        mol._x = (mol._x + random.uniform(-b, b)) % self._lx
        mol._y = (mol._y + random.uniform(-b, b)) % self._ly
        mol._z = (mol._z + random.uniform(-b, b)) % self._lz

        # Calculate energy after the move
        new_energy = self._energy_of_molecule(i)

        delta_E = new_energy - old_energy

        # acept move if ΔE <= 0, else with probability exp(-delta_E / T)
        if delta_E <= 0 or random.random() < math.exp(-delta_E / T):
            return True  # Move accepted
        else:
            # Restore original position (move rejected)
            mol._x, mol._y, mol._z = old_x, old_y, old_z
            return False

    def _energy_of_molecule(self, idx):
        """
        Calculate the energy of one molecule (idx) with all other molecules in the box..
        """
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

            # Apply periodic boundary conditions (minimum image)
            dx -= round(dx / self._lx) * self._lx
            dy -= round(dy / self._ly) * self._ly
            dz -= round(dz / self._lz) * self._lz

            r_squared = dx * dx + dy * dy + dz * dz

            if r_squared < cutoff_squared and r_squared > 1e-12:
                inv_r2 = 1.0 / r_squared
                inv_r6 = inv_r2 ** 3
                inv_r12 = inv_r6 ** 2
                energy += 4 * (inv_r12 - inv_r6) + shift_value
        return energy
