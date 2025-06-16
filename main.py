from src.ljts.boks import Box
import matplotlib.pyplot as plt

def main():
    length_x = 5
    length_y = 40
    length_z = 5
    density_liquid = 0.73
    density_vapor = 0.02

    box = Box(length_x, length_y, length_z, density_liquid, density_vapor)

    initial_energy = box.compute_potential_energy()
    print(f"Potential energy FØR MC: {initial_energy:.6f}")

    N = box.get_molecule_count()
    b = 1.0 / 8
    T = 0.8  
    sweeps = 50
    log_interval = 1

    energy_history = [initial_energy]
    accept_total = 0

    for sweep in range(1, sweeps+1):
        accept = box.mc_sweep(b, T)
        accept_total += accept

        if sweep % log_interval == 0:
            curr_energy = box.compute_potential_energy()
            energy_history.append(curr_energy) #maybe use for plotting later if i have time
            print(f"Sweep {sweep}: Pot.energy = {curr_energy:.6f}")

    accept_rate = accept_total / (sweeps * N)
    print(f"\nAcceptance ratio: {accept_rate:.3f}")


if __name__ == "__main__":
    main()

