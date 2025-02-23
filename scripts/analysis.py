import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.fft import fft, fftfreq

# Define paths
base_path = r"C:\Users\monica\Desktop\lammps_tutorial\polymer_cyclic_loading\results"
output_path = base_path  # Save outputs in the results folder

# Time conversion factor (LAMMPS timestep = 0.001 LJ units → convert to ns)
timestep_size = 0.001  # Assuming LAMMPS uses 0.001 LJ units
time_factor = timestep_size * 1e3  # Convert to nanoseconds (ns)

# Helper function to safely load data, ignoring comments
def safe_loadtxt(file_path, usecols=None, skip_metadata_rows=0):
    try:
        data = []
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue  # Skip comments
                values = line.split()
                if len(values) < max(usecols) + 1:
                    continue  # Skip malformed rows
                data.append([float(values[i]) for i in usecols])
        return np.array(data)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

# 1. **Plot End-to-End Distance Distribution**
def plot_end_to_end():
    file_path = os.path.join(base_path, "end_to_end_distribution.txt")
    data = safe_loadtxt(file_path, usecols=[1, 2], skip_metadata_rows=3)  # Bin position, count

    if data is None or len(data) == 0:
        print("No valid data for End-to-End distance.")
        return

    bins, values = data[:, 0], data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.bar(bins, values, width=bins[1] - bins[0], alpha=0.7, color='b', label="End-to-End Distance")
    plt.xlabel("End-to-End Distance")
    plt.ylabel("Frequency")
    plt.title("End-to-End Distance Distribution")
    plt.legend()
    plt.grid()

    output_file = os.path.join(output_path, "end_to_end_distribution.png")
    plt.savefig(output_file)
    plt.show()
    print(f"Saved: {output_file}")

# 2. **Plot Radius of Gyration (Rg) Evolution**
def plot_radius_of_gyration():
    file_path = os.path.join(base_path, "radius_of_gyration.txt")

    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
        
        # Extracting valid data (skipping metadata rows)
        time_steps = []
        rg_values = []
        i = 0
        while i < len(lines):
            if lines[i].startswith("#") or not lines[i].strip():
                i += 1
                continue  # Skip headers or empty lines
            
            # Ensure we have a valid pair of lines
            if i + 1 < len(lines):
                time_step = float(lines[i].split()[0])  # Extract timestep
                rg_value = float(lines[i + 1].split()[1])  # Extract Rg value
                time_steps.append(time_step)
                rg_values.append(rg_value)
                i += 2  # Move to the next timestep
            else:
                break

        # Convert time to nanoseconds
        time_steps = np.array(time_steps) * time_factor
        rg_values = np.array(rg_values)

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    plt.figure(figsize=(8, 6))
    plt.plot(time_steps, rg_values, marker='o', linestyle='-', color='r', label="Radius of Gyration")
    plt.xlabel("Time (ns)")
    plt.ylabel("Radius of Gyration (Rg)")
    plt.title("Evolution of Radius of Gyration")
    plt.legend()
    plt.grid()

    output_file = os.path.join(output_path, "radius_of_gyration.png")
    plt.savefig(output_file)
    plt.show()
    print(f"Saved: {output_file}")

# 3. **Compute and Plot Stress Autocorrelation Function**
def plot_stress_and_moduli():
    file_path = os.path.join(base_path, "stress_correlation.txt")

    try:
        with open(file_path, "r") as f:
            lines = f.readlines()

        # Locate the actual data section
        data_start_idx = None
        for i, line in enumerate(lines):
            if not line.startswith("#") and line.strip():
                data_start_idx = i
                break

        if data_start_idx is None:
            raise ValueError("No valid data found in stress_correlation.txt")

        # Load data skipping headers
        data = []
        for line in lines[data_start_idx:]:
            parts = line.split()
            if len(parts) >= 4 and parts[1].isdigit():  # Ensure proper numeric values
                time_delta = float(parts[1])  # Second column is TimeDelta
                autocorr = float(parts[3])  # Fourth column is v_pxy*v_pxy
                data.append((time_delta, autocorr))

        data = np.array(data)
        if data.size == 0:
            raise ValueError("No valid stress autocorrelation data found")

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    # Extract time and stress autocorrelation
    time_steps, stress_values = data[:, 0] * time_factor, data[:, 1]  # Convert to ns

    # Plot Stress Autocorrelation
    plt.figure(figsize=(8, 6))
    plt.plot(time_steps, stress_values, marker='s', linestyle='-', color='g', label="Stress Autocorrelation")
    plt.xlabel("Time (ns)")
    plt.ylabel("Stress Autocorrelation")
    plt.title("Stress Autocorrelation Function")
    plt.legend()
    plt.grid()

    output_file = os.path.join(output_path, "stress_autocorrelation.png")
    plt.savefig(output_file)
    plt.show()
    print(f"Saved: {output_file}")

    # ---- Compute Storage Modulus (G') and Loss Modulus (G'') ----
    dt = np.mean(np.diff(time_steps))  # Time step size (assumed uniform)
    omega = np.fft.rfftfreq(len(time_steps), d=dt) * 2 * np.pi  # Angular frequencies

    fft_result = np.fft.rfft(stress_values) * dt  # Fourier transform of stress autocorrelation
    G_prime = omega * np.real(fft_result)  # Storage modulus G'
    G_double_prime = omega * np.imag(fft_result)  # Loss modulus G''

    # ---- Plot G' and G'' ----
    plt.figure(figsize=(8, 6))
    plt.loglog(omega[1:], G_prime[1:], marker='o', linestyle='-', label="G' (Storage Modulus)", color='b')
    plt.loglog(omega[1:], G_double_prime[1:], marker='s', linestyle='-', label="G'' (Loss Modulus)", color='r')

    plt.xlabel("Frequency ω (rad/ns)")
    plt.ylabel("Modulus (G', G'') ")
    plt.title("Storage and Loss Moduli")
    plt.legend()
    plt.grid()

    output_file = os.path.join(output_path, "moduli_Gp_Gpp.png")
    plt.savefig(output_file)
    plt.show()
    print(f"Saved: {output_file}")

# Run all analysis functions
if __name__ == "__main__":
    try:
        plot_end_to_end()
    except Exception as e:
        print(f"Error plotting End-to-End Distance: {e}")

    try:
        plot_radius_of_gyration()
    except Exception as e:
        print(f"Error plotting Radius of Gyration: {e}")

    try:
        plot_stress_and_moduli()
    except Exception as e:
        print(f"Error plotting stress and moduli: {e}")

    print("Analysis complete. Check the generated PNG files in the results folder.")
