"""
Multi-step pathway visualization with arrows and reagents.
- Each stage is a molecule (SMILES).
- Arrows are drawn between molecules.
- Reagents/conditions are displayed above arrows.
"""

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os

# -------------------------------
# Setup
# -------------------------------
os.makedirs("data", exist_ok=True)

# -------------------------------
# 1. Define stages and reagents
# -------------------------------
stage_smiles = [
    "C1(=CC=C(C=C1C=CC(=O)OC)OC)OC",  # Stage 1
    "C(C)OCC",                         # Stage 2 (example intermediate)
    "C1(=CC=C(C=C1C2C(C2)C(=O)OC)OC)OCC1=CC(=CC=C1CC(=O)O)Cl",  # Stage 3
    "C1(=CC=CC=C1)C",                  # Stage 4
    "C1=CC(=CC=C1CC(=O)Cl)Cl"          # Final product
]

# Reagents above arrows (between stages)
reagents = [
    "Pd(II), ClCCl",
    "Base, Heat",
    "Amine, DMF",
    "SOCl2"
]

# -------------------------------
# 2. Convert SMILES to molecules
# -------------------------------
mols = [Chem.MolFromSmiles(s) for s in stage_smiles]

# Compute 2D coordinates
for mol in mols:
    Chem.rdDepictor.Compute2DCoords(mol)

# -------------------------------
# 3. Draw each molecule as image
# -------------------------------
subImgSize = (250, 250)
mol_imgs = [Draw.MolToImage(mol, size=subImgSize) for mol in mols]

# -------------------------------
# 4. Combine into one canvas with arrows + labels
# -------------------------------
arrow_width = 120
canvas_width = len(mols) * subImgSize[0] + (len(mols)-1) * arrow_width
canvas_height = subImgSize[1] + 80  # extra space for reagent labels

canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
draw = ImageDraw.Draw(canvas)

# Font for reagent labels
try:
    font = ImageFont.truetype("arial.ttf", 20)
except:
    font = ImageFont.load_default()

# Paste molecules and draw arrows with labels
x = 0
for i, mol_img in enumerate(mol_imgs):
    # Paste molecule
    canvas.paste(mol_img, (x, 50))
    x += subImgSize[0]

    # Draw arrow + reagent label if not the last molecule
    if i < len(reagents):
        arrow_start = x
        arrow_end = x + arrow_width - 20
        y = subImgSize[1] // 2 + 50

        # Draw arrow line
        draw.line((arrow_start, y, arrow_end, y), fill="black", width=3)
        # Draw arrowhead
        draw.polygon([(arrow_end, y), (arrow_end-10, y-5), (arrow_end-10, y+5)], fill="black")

        # Draw reagent label above arrow
        bbox = draw.textbbox((0, 0), reagents[i], font=font)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]
        text_x = (arrow_start + arrow_end) // 2 - text_w // 2
        text_y = y - text_h - 10  # 10 px above arrow
        draw.text((text_x, text_y), reagents[i], font=font, fill="black")

        # Advance x for next molecule
        x += arrow_width

# -------------------------------
# 5. Save final pathway image
# -------------------------------
canvas.save("data/multi_step_pathway.png")
print("Pathway saved as data/multi_step_pathway.png")
