from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os

"""
Draws a multi-step synthesis pathway from a list of SMIRKS-like strings:
Format of each step: "reactants > reagents > products"
Example:
[
    "CCO>H2SO4>CCOC",
    "CCOC>NaOH>CC=O"
]
"""

def center_text(draw, text, x_center, y, font):
    """Draw text centered horizontally."""
    bbox = draw.textbbox((0, 0), text, font=font)
    w = bbox[2] - bbox[0]
    h = bbox[3] - bbox[1]
    draw.text((x_center - w // 2, y), text, font=font, fill="black")


def draw_multistep_pathway(smirk_steps, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Parse steps into reactants, reagents, products
    pathway = []
    for step in smirk_steps:
        parts = step.split(">")
        if len(parts) != 3:
            raise ValueError(f"Invalid step (must be reactants>reagents>products): {step}")

        reactants = parts[0].split(".") if parts[0] else []
        reagents = parts[1].split(".") if parts[1] else []
        products = parts[2].split(".") if parts[2] else []
        pathway.append((reactants, reagents, products))

    # Gather all molecule SMILES in drawing order
    all_smiles = []
    for reactants, _, products in pathway:
        all_smiles.extend(reactants)
        all_smiles.extend(products)

    # Convert to RDKit mols
    mols = []
    for sm in all_smiles:
        mol = Chem.MolFromSmiles(sm)
        if mol:
            Chem.rdDepictor.Compute2DCoords(mol)
            mols.append(mol)
        else:
            print(f" Warning: Could not parse SMILES: {sm}")

    subImgSize = (130, 130)
    mol_imgs = [Draw.MolToImage(m, size=subImgSize) for m in mols]

    arrow_width = 140
    canvas_width = len(mol_imgs) * (subImgSize[0] + 30) + (len(pathway) - 1) * arrow_width
    canvas_height = subImgSize[1] + 140

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    try:
        font = ImageFont.truetype("arial.ttf", 18)
    except:
        font = ImageFont.load_default()

    x = 0
    img_index = 0

    for step_i, (reactants, reagents, products) in enumerate(pathway):
        num_molecules = len(reactants) + len(products)

        for mol_i in range(num_molecules):
            # Draw molecule
            canvas.paste(mol_imgs[img_index], (x, 60))

            # Label molecule (SMILES for now — can replace w/IUPAC)
            mol_name = Chem.MolToSmiles(mols[img_index], canonical=True)
            center_text(draw, mol_name, x + subImgSize[0] // 2, 60 + subImgSize[1] + 8, font)

            # Draw "+" between molecules in the same step
            if mol_i < num_molecules - 1:
                plus_x = x + subImgSize[0] + 10
                plus_y = 60 + subImgSize[1] // 2 - 10
                draw.text((plus_x, plus_y), "+", font=font, fill="black")
                x += subImgSize[0] + 40
            else:
                x += subImgSize[0] + 10

            img_index += 1

        # Draw arrow & reagent label if not last step
        if step_i < len(pathway) - 1:
            y = subImgSize[1] // 2 + 60
            arrow_start = x
            arrow_end = x + arrow_width - 20

            # Arrow line + head
            draw.line((arrow_start, y, arrow_end, y), fill="black", width=4)
            draw.polygon([(arrow_end, y), (arrow_end - 12, y - 6), (arrow_end - 12, y + 6)], fill="black")

            # Center reagent label above arrow
            arrow_center = (arrow_start + arrow_end) // 2
            center_text(draw, ", ".join(reagents), arrow_center, y - 45, font)

            x += arrow_width

    canvas.save(output_path)
    print(f" Multistep synthesis pathway saved to: {output_path}")


# Example usage
if __name__ == "__main__":
    steps = [
        "CCO>H2SO4>CCOC",
        "CCOC>NaOH>CC=O",
        "CC=O>KMnO4>CC(=O)O"
    ]
    
    draw_multistep_pathway(
        steps,
        os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "data", "Script_4", "attempted_synthesis.png"))
    )
   #draw_multistep_pathway(steps, "src/data/attempted_synthesis.png")
