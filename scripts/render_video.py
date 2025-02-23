from ovito.io import import_file
from ovito.vis import Viewport
from ovito.pipeline import Pipeline
from ovito.modifiers import ColorCodingModifier, CreateBondsModifier
import sys
import os

# Get input file from command-line argument
if len(sys.argv) < 2:
    print("Usage: py render_video.py <path_to_lammpstrj>")
    sys.exit(1)

input_file = sys.argv[1]

# Ensure output is saved in the same folder as the LAMMPS trajectory
output_folder = os.path.dirname(input_file)
output_file = os.path.join(output_folder, os.path.basename(input_file).replace(".lammpstrj", ".mp4"))

# Load trajectory
pipeline = import_file(input_file)
pipeline.add_to_scene()

# Ensure atoms are visible by adding color coding (optional)
pipeline.modifiers.append(ColorCodingModifier(property='Particle Type'))

# If bonds are missing, generate them (adjust cutoff as needed)
bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
bonds_modifier.cutoff = 1.5  # Adjust cutoff based on your polymer model
pipeline.modifiers.append(bonds_modifier)

# Configure viewport
vp = Viewport()
vp.type = Viewport.Type.PERSPECTIVE  # Options: PERSPECTIVE, TOP, FRONT, SIDE
vp.camera_pos = (0, 0, -50)  # Adjust camera position
vp.camera_dir = (0, 0, 1)  # Ensure camera is pointing in the right direction

# Render animation with a black background
vp.render_anim(size=(1920, 1080), filename=output_file, fps=30, background=(0, 0, 0))

print(f"Video saved as {output_file}")
