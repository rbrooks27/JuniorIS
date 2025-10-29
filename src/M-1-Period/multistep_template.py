from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os

'''
Overview of this script:
This script provides a function to draw multistep chemical reaction pathways from a single SMIRKS
line formatted as Reactants > Reagents > Products. It uses RDKit to parse the SMILES and PIL to create
a visual representation of the reaction flow, including arrows and reagent labels.
'''

def draw_multistep_pathway_from_smirks(smirks, output_path):
    """
    Draws a multistep pathway from a single SMIRKS line.
    Format: reactants > reagents > products
    """
    os.makedirs("data", exist_ok=True)

    # Parse SMIRKS line
    parts = smirks.split(">")
    if len(parts) != 3:
        raise ValueError(f"Invalid SMIRKS format: {smirks}")

    reactants = parts[0].split(".") if parts[0] else []
    reagents = parts[1].split(".") if parts[1] else []
    products = parts[2].split(".") if parts[2] else []

    # Define stages (reactants + products)
    stage_smiles = reactants + products

    # Convert SMILES to molecules
    mols = [Chem.MolFromSmiles(s) for s in stage_smiles if s]
    for mol in mols:
        if mol:
            Chem.rdDepictor.Compute2DCoords(mol)

    # Draw molecules
    subImgSize = (100, 100)
    mol_imgs = [Draw.MolToImage(mol, size=subImgSize) for mol in mols if mol]

    # Canvas size
    arrow_width = 120
    canvas_width = len(mol_imgs) * subImgSize[0] + (len(mol_imgs)-1) * arrow_width
    canvas_height = subImgSize[1] + 80

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Font
    try:
        font = ImageFont.truetype("arial.ttf", 20)
    except:
        font = ImageFont.load_default()

    # Paste molecules + arrows
    x = 0
    for i, mol_img in enumerate(mol_imgs):
        canvas.paste(mol_img, (x, 50))
        x += subImgSize[0]

        # Draw arrow with reagents above (only once, between reactants and products)
        if i == len(reactants) - 1 and reagents:
            arrow_start = x
            arrow_end = x + arrow_width - 20
            y = subImgSize[1] // 2 + 50
            draw.line((arrow_start, y, arrow_end, y), fill="black", width=3)
            draw.polygon([(arrow_end, y), (arrow_end-10, y-5), (arrow_end-10, y+5)], fill="black")

            # Join reagents into label
            reagents_label = ", ".join(reagents)
            bbox = draw.textbbox((0, 0), reagents_label, font=font)
            text_w = bbox[2] - bbox[0]
            text_h = bbox[3] - bbox[1]
            text_x = (arrow_start + arrow_end) // 2 - text_w // 2
            text_y = y - text_h - 10
            draw.text((text_x, text_y), reagents_label, font=font, fill="black")

            x += arrow_width

    # Save
    canvas.save(output_path)
    print(f"Pathway saved as {output_path}")
