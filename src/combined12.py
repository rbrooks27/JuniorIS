from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os

# -------------------------------
# Helper functions
# -------------------------------
def parse_smirks_line(smirks_line):
    """
    Parse a SMIRKS line into a list of (reactant, reagent, product) tuples.
    Handles multiple reactants/products and canonicalizes SMILES.
    Picks main product by largest heavy atom count.
    """
    parts = smirks_line.strip().split(">")
    if len(parts) != 3:
        raise ValueError(f"Invalid SMIRKS format: {smirks_line}")

    reactants = [Chem.MolToSmiles(Chem.MolFromSmiles(s), canonical=True) for s in parts[0].split(".") if s]
    reagents = [s.strip() for s in parts[1].split(".") if s]
    products = [Chem.MolToSmiles(Chem.MolFromSmiles(s), canonical=True) for s in parts[2].split(".") if s]

    # Pick main product (largest heavy atom count)
    def main_product(smiles_list):
        return max(smiles_list, key=lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms())

    main_prod = main_product(products)

    # Combine reactants with main product into step tuples
    steps = []
    for react in reactants:
        steps.append((react, ", ".join(reagents), main_prod))

    return steps

# -------------------------------
# Drawing function
# -------------------------------
def draw_multistep_pathway(steps, output_path):
    """
    Draws a clean multistep pathway with molecules, arrows, and reagents.
    `steps`: list of tuples (reactant, reagent, product)
    """
    if not steps:
        raise ValueError("No steps to draw.")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Extract stage molecules and reagents
    stage_smiles = [steps[0][0]] + [p for (_, _, p) in steps]
    reagents = [r for (_, r, _) in steps]

    # Convert to RDKit molecules and compute 2D coordinates
    mols = [Chem.MolFromSmiles(s) for s in stage_smiles]
    for mol in mols:
        Chem.rdDepictor.Compute2DCoords(mol)

    # Generate molecule images
    subImgSize = (250, 250)
    mol_imgs = [Draw.MolToImage(mol, size=subImgSize) for mol in mols]

    # Canvas dimensions
    arrow_width = 120
    canvas_width = len(mols) * subImgSize[0] + (len(mols)-1) * arrow_width
    canvas_height = subImgSize[1] + 80

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Font for reagent labels
    try:
        font = ImageFont.truetype("arial.ttf", 20)
    except:
        font = ImageFont.load_default()

    # Paste molecules and draw arrows
    x = 0
    for i, mol_img in enumerate(mol_imgs):
        canvas.paste(mol_img, (x, 50))
        x += subImgSize[0]

        if i < len(reagents):
            # Arrow
            arrow_start = x
            arrow_end = x + arrow_width - 20
            y = subImgSize[1] // 2 + 50
            draw.line((arrow_start, y, arrow_end, y), fill="black", width=3)
            draw.polygon([(arrow_end, y), (arrow_end-10, y-5), (arrow_end-10, y+5)], fill="black")

            # Reagent label
            bbox = draw.textbbox((0, 0), reagents[i], font=font)
            text_w = bbox[2] - bbox[0]
            text_h = bbox[3] - bbox[1]
            text_x = (arrow_start + arrow_end) // 2 - text_w // 2
            text_y = y - text_h - 10
            draw.text((text_x, text_y), reagents[i], font=font, fill="black")

            x += arrow_width

    # Save the image
    canvas.save(output_path)
    print(f"Pathway saved as {output_path}")


# -------------------------------
# Example usage
# -------------------------------
if __name__ == "__main__":
    smirks_example = "CCO>H2SO4>C=C.CCO>Br2>Br-CH-CH3.Br-CH-CH3>NH3>CH3CH2NH2" # replace with any SMIRKS line from your dataset
    steps = parse_smirks_line(smirks_example)
    output_file = os.path.join("data", "multi_step_pathway.png")
    draw_multistep_pathway(steps, output_file)
