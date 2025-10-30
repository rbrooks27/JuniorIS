import os
from tqdm import tqdm
from combined123 import parse_multistep_smirks, draw_multistep_pathway  #  unified script
from reactionGraph import ReactionGraph # this for graphing structures

"""
Overview:
This script reads a dataset of chemical reactions in SMIRKS format (Reactants > Reagents > Products)
and generates visual representations of the reaction pathways using the unified parser/drawing functions.
"""

# -------------------------------
# Dataset and output settings
# -------------------------------
dataset_file = "dataset/textbookReactions.txt"  # smaller sample for testing
# dataset_file = "dataset/reactionSmilesFigShare2024.txt"  # full dataset

output_dir = "data/reactions"  # output directory for pathway images
os.makedirs(output_dir, exist_ok=True)

MAX_LINES = 10  # set to None to process all lines

# Initialization of Reaction Graph
G = ReactionGraph()

# -------------------------------
# Process each SMIRKS line
# -------------------------------
with open(dataset_file, "r") as f:
    for i, line in enumerate(tqdm(f, desc="Generating pathways and building graph")):
        if MAX_LINES and i >= MAX_LINES:
            break

        line = line.strip()
        if not line or ">" not in line:
            continue

        try:
            # Parse the SMIRKS into steps
            steps = parse_multistep_smirks(line)
            # Skip if parsing returned nothing
            if not steps:
                print(f"No valid steps found at line {i}")
                continue
            
            G.load_steps(steps)  # Load steps into the reaction graph

            # Generate output path
            output_path = os.path.join(output_dir, f"reaction_{i}.png")

            # Draw the multistep pathway
            draw_multistep_pathway(steps, output_path)

        except Exception as e:
            print(f"Error processing line {i}: {line}")
            print(f"Error: {e}")
            continue

print("Done generating reaction pathways!")

starting_materials = {"CCO", "H2SO4", "C=C"} # Example starting materials
target = "OCC" # Example target molecule

# Find a synthesis path in the reaction graph
path = G.find_path(starting_materials, target)

if path:
    print("Synthesis path found:")
    for react, reagents, prod in path:
        print(f"{react} --[{reagents}]--> {prod}")
    
    final_path_img = os.path.join(output_dir, "synthesis_pathway.png")
    draw_multistep_pathway(path, final_path_img)
    print(f"Synthesis pathway image saved as {final_path_img}")
else:
    print("No synthesis path found from starting materials to target.")