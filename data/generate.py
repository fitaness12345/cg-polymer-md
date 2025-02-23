import math, random

# User parameters
num_atoms = 1000
box_size = 30.0
target_bond_length = 0.96  # Approximate equilibrium bond length
max_attempts = 1000

# Start at a random position well inside the box
margin = 5.0
x = random.uniform(margin, box_size - margin)
y = random.uniform(margin, box_size - margin)
z = random.uniform(margin, box_size - margin)
positions = [(x, y, z)]

# Generate a random walk with self-avoidance
for i in range(1, num_atoms):
    attempt = 0
    while attempt < max_attempts:
        attempt += 1
        # Random direction in 3D: sample from spherical coordinates
        theta = random.uniform(0, math.pi)
        phi = random.uniform(0, 2*math.pi)
        dx = target_bond_length * math.sin(theta) * math.cos(phi)
        dy = target_bond_length * math.sin(theta) * math.sin(phi)
        dz = target_bond_length * math.cos(theta)
        new_x = positions[-1][0] + dx
        new_y = positions[-1][1] + dy
        new_z = positions[-1][2] + dz
        # Check if new position is within the box
        if new_x < 0 or new_x > box_size or new_y < 0 or new_y > box_size or new_z < 0 or new_z > box_size:
            continue
        # Check against previous positions for overlap
        too_close = False
        for (px, py, pz) in positions:
            dist_sq = (new_x-px)**2 + (new_y-py)**2 + (new_z-pz)**2
            if dist_sq < (0.8 * target_bond_length)**2:  # some tolerance
                too_close = True
                break
        if not too_close:
            positions.append((new_x, new_y, new_z))
            break
    if attempt == max_attempts:
        print("Failed to place atom", i)
        break

# Write LAMMPS data file
with open("initial_polymer_data.dat", "w") as f:
    f.write("LAMMPS Description\n\n")
    f.write(f"{len(positions)} atoms\n")
    f.write(f"{len(positions)-1} bonds\n\n")
    f.write("1 atom types\n")
    f.write("1 bond types\n\n")
    f.write(f"0.0 {box_size} xlo xhi\n")
    f.write(f"0.0 {box_size} ylo yhi\n")
    f.write(f"0.0 {box_size} zlo zhi\n\n")
    f.write("Masses\n\n")
    f.write("1 1.0\n\n")
    f.write("Atoms\n\n")
    for i, (x, y, z) in enumerate(positions):
        # all atoms belong to molecule 1
        f.write(f"{i+1} 1 1 {x:.3f} {y:.3f} {z:.3f}\n")
    f.write("\nBonds\n\n")
    for i in range(1, len(positions)):
        f.write(f"{i} 1 {i} {i+1}\n")

print("Initial polymer data generated in 'initial_polymer_data.dat'")
