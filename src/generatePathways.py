# just general imports needed
import os
from tqdm import tqdm
from rdkit import Chem

# Functions from combined123.py and reactionGraph.py which are used for graph building and pathway drawing
from combined123 import parse_multistep_smirks, draw_multistep_pathway  #  unified script
from reactionGraph import ReactionGraph


"""
# Overview:
# This script reads a dataset of chemical reactions in SMIRKS format (Reactants > Reagents > Products)
# and generates visual representations of the reaction pathways using the unified parser/drawing functions.
"""

# --------------------------------------------
# Helpers
# --------------------------------------------
def canon(smiles: str) -> str:
    """Convert SMILES to canonical form."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return smiles


# --------------------------------------------
# Dataset + Output
# --------------------------------------------
dataset_file = "dataset/textbookReactions.txt"   # small dataset for testing
# dataset_file = "dataset/reactionSmilesFigShare2024.txt"  # full dataset

output_dir = "data/reactions"
os.makedirs(output_dir, exist_ok=True)

MAX_LINES = 20  # set None to process entire file

# Initialize Reaction Graph
G = ReactionGraph()


# --------------------------------------------
# Build graph + generate images
# --------------------------------------------
with open(dataset_file, "r") as f:
    for i, line in enumerate(tqdm(f, desc="Generating pathways and building graph")):
        if MAX_LINES and i >= MAX_LINES:
            break

        line = line.strip()
        if not line or ">" not in line:
            continue

        try:
            # Parse the SMIRKS into individual (reactant, reagent, product) steps
            steps = parse_multistep_smirks(line)

            if not steps:
                print(f"No valid steps at line {i}: {line}")
                continue

            # Store steps in the reaction graph
            G.load_steps(steps)

            # Output pathway drawing
            output_path = os.path.join(output_dir, f"reaction_{i}.png")
            draw_multistep_pathway(steps, output_path)

        except Exception as e:
            print(f"Error processing line {i}: {line}")
            print(f"Error: {e}")
            continue

print("Done generating individual reaction pathways.")


# --------------------------------------------
# Retrosynthesis Search
# --------------------------------------------
# Canonicalize all starting materials
# starting_materials_raw = ["CCO"]
# starting_materials_raw = ["CC=O"]
starting_materials_raw = ["C1(=CC=C(C=C1C[P](=O)(OCC)OCC)[S](F)(F)(F)(F)F)[N+](=O)[O-]", "C1(=CC=CC=C1)C=O>N(C)(C)C=O.CC(C)(C)[O-]", "[K+]"]
starting_materials = {canon(s) for s in starting_materials_raw}

# Canonicalize target!
# target_raw = "COCC"
# target_raw = "CC(=O)O"
target_raw = "C1(=CC=C(C=C1C=CC2=CC=CC=C2)[S](F)(F)(F)(F)F)[N+](=O)[O-]"
target = canon(target_raw)

print(f"\nSearching for path to target: {target_raw} (canonical: {target})")

path = G.find_path(starting_materials, target)


# --------------------------------------------
# Output Retrosynthesis Path
# --------------------------------------------
if path:
    print("\nSynthesis path found:\n")
    for react, reagents, prod in path:
        print(f"{react} --[{reagents}]--> {prod}")

    final_path_img = os.path.join(output_dir, "synthesis_pathway.png")
    draw_multistep_pathway(path, final_path_img)
    print(f"\nSynthesis pathway image saved as:\n{final_path_img}\n")

else:
    print("\nNo synthesis path found from starting materials to target.\n")

