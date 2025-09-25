import os
from multistep_template import draw_multistep_pathway_from_smirks
from tqdm import tqdm  # progress bar
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

dataset_file = "dataset/reactionSmilesFigShare2024.txt"   # large dataset
output_dir = "data/reactions"
os.makedirs(output_dir, exist_ok=True)

# Optional: limit processing for testing
MAX_LINES = 10  # set to 1000 for testing, or None for all

# Read dataset and generate pathways
with open(dataset_file, "r") as f:
    # Progress bar for large files
    for i, line in enumerate(tqdm(f, desc="Generating pathways")):
        # Optional line limit for testing
        if MAX_LINES and i >= MAX_LINES:
            break

        line = line.strip()
        if not line:
            continue

        try:
            # Parse line directly as SMIRKS (Reactants > Reagents > Products)
            rxn = AllChem.ReactionFromSmarts(line)
            if rxn is None:
                print(f"Skipping invalid reaction at line {i}")
                continue

            # Draw reaction image
            img = Draw.ReactionToImage(rxn, subImgSize=(250, 250))
            output_path = os.path.join(output_dir, f"reaction_{i}.png")
            img.save(output_path)
        except Exception as e:
            print(f"Skipping badly formatted line {i}: {line}")
            continue



'''
for filename in os.listdir(dataset_dir):
    if not filename.endswith(".txt"):  # only parse text dataset files
        continue

    filepath = os.path.join(dataset_dir, filename)
    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            
            # Example line format: SMILES1,SMILES2,...; Reagent1,Reagent2,...
            try:
                stages_part, reagents_part = line.split(";")
                stage_smiles = [s.strip() for s in stages_part.split(",")]
                reagents = [r.strip() for r in reagents_part.split(",")]
            except ValueError:
                print(f"Skipping badly formatted line in {filename}: {line}")
                continue

            output_path = os.path.join(output_dir, f"{filename}_pathway_{i}.png")
            draw_multistep_pathway(stage_smiles, reagents, output_path)
'''