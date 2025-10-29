from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os

'''
Overview of this script:
This script provides a function to draw multistep chemical reaction pathways
from a single or multistep SMIRKS string formatted as:
Reactant > Reagent > Product; Product > Reagent > Product2; ...
It uses RDKit to parse SMILES and PIL to visualize the reaction flow.
'''

def smirks_to_stages(smirks):
    """
    Convert a multi-step SMIRKS string into stage_smiles and reagents lists.

    Example:
        Input:  "CCO>H2SO4>CC=O; CC=O>NaBH4>CCO"
        Output: stage_smiles = ["CCO", "CC=O", "CCO"]
                reagents = ["H2SO4", "NaBH4"]
    """
    steps = [step.strip() for step in smirks.split(";") if step.strip()]
    stage_smiles = []
    reagents = []

    for i, step in enumerate(steps):
        try:
            reactant, reagent, product = [x.strip() for x in step.split(">")]
        except ValueError:
            print(f"[Warning] Invalid SMIRKS skipped: {step}")
            continue

        if i == 0:
            stage_smiles.append(reactant)

        stage_smiles.append(product)
        reagents.append(reagent)

    # Validate with RDKit (remove invalid SMILES)
    valid_stage_smiles = []
    for s in stage_smiles:
        mol = Chem.MolFromSmiles(s)
        if mol:
            valid_stage_smiles.append(s)
        else:
            print(f"[Warning] Invalid SMILES skipped: {s}")

    return valid_stage_smiles, reagents


def draw_multistep_pathway_from_smirks(smirks, output_path="data/msp.png"):
    """
    Draws a multistep pathway from a SMIRKS string.
    Example:
        "CCO>H2SO4>CC=O; CC=O>NaBH4>CCO"
    """
    os.makedirs("data", exist_ok=True)

    # Parse into stages
    stage_smiles, reagents = smirks_to_stages(smirks)
    if not stage_smiles:
        raise ValueError("No valid SMILES structures found in input SMIRKS.")

    # Convert SMILES to RDKit molecules
    mols = [Chem.MolFromSmiles(s) for s in stage_smiles]
    for mol in mols:
        Chem.rdDepictor.Compute2DCoords(mol)

    # Generate molecule images
    subImgSize = (150, 150)
    mol_imgs = [Draw.MolToImage(mol, size=subImgSize) for mol in mols]

    # Canvas layout
    arrow_width = 140
    canvas_width = len(mol_imgs) * subImgSize[0] + (len(mol_imgs) - 1) * arrow_width
    canvas_height = subImgSize[1] + 100

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Font setup
    try:
        font = ImageFont.truetype("arial.ttf", 18)
    except:
        font = ImageFont.load_default()

    # Draw molecules and arrows
    x = 0
    for i, mol_img in enumerate(mol_imgs):
        canvas.paste(mol_img, (x, 60))
        x += subImgSize[0]

        # Draw arrow and reagent between stages
        if i < len(reagents):
            arrow_start = x
            arrow_end = x + arrow_width - 30
            y = subImgSize[1] // 2 + 60
            draw.line((arrow_start, y, arrow_end, y), fill="black", width=3)
            draw.polygon([(arrow_end, y), (arrow_end - 10, y - 5), (arrow_end - 10, y + 5)], fill="black")

            reagent_label = reagents[i]
            bbox = draw.textbbox((0, 0), reagent_label, font=font)
            text_w = bbox[2] - bbox[0]
            text_h = bbox[3] - bbox[1]
            text_x = (arrow_start + arrow_end) // 2 - text_w // 2
            text_y = y - text_h - 15
            draw.text((text_x, text_y), reagent_label, font=font, fill="black")

            x += arrow_width

    # Save output image
    canvas.save(output_path)
    print(f"Pathway saved to {output_path}")


# -----------------------------
# Example usage
# -----------------------------
if __name__ == "__main__":
    smirks = "CCO>H2SO4>CC=O; CC=O>NaBH4>CCO; CCO>H2O>CC(=O)O"
    draw_multistep_pathway_from_smirks(smirks)
