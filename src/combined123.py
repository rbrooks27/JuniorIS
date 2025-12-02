from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os

# -------------------------------
# Helper functions
# -------------------------------
def parse_multistep_smirks(smirks_string, separator=";"):
    """
    Parse a string containing multiple reactions separated by `separator` (default ';' or '.').
    Returns a list of (reactant, reagent, product) tuples.
    """
    # Clean up input
    smirks_string = smirks_string.strip()
    
    # Split the string into individual reactions
    # Allow both ';' or '.' as separators
    # raw_reactions = [s.strip() for s in smirks_string.replace(".", separator).split(separator) if ">" in s]
    
    raw_reactions = [s.strip() for s in smirks_string.split(separator) if ">" in s]
    
    all_steps = []

    for rxn_line in raw_reactions:
        parts = rxn_line.split(">")
        if len(parts) != 3:
            continue # skip malformed
            # raise ValueError(f"Invalid SMIRKS format: {rxn_line}")
        reactant_block, reagent_block, product_block = parts

        # Old reactants, reagents, products parsing
        """
        reactants = [Chem.MolToSmiles(Chem.MolFromSmiles(s), canonical=True) for s in parts[0].split(".") if s]
        reagents = [s.strip() for s in parts[1].split(".") if s]
        products = [Chem.MolToSmiles(Chem.MolFromSmiles(s), canonical=True) for s in parts[2].split(".") if s]
        """
        # Reactants: split by '.'
        reactants = [s for s in reactant_block.split(".") if s]
        
        # Reagents: keep as a single string (commas within reagents allowed)
        # reagents = [s for s in reagent_block.split(".") if s]
        reagents = ",".join([s for s in reagent_block.split(".") if s])
        
        products = []
        for s in product_block.split("."):
            if not s:
                continue
            mol = Chem.MolFromSmiles(s)
            if mol:
                products.append(Chem.MolToSmiles(mol, canonical=True))
            else:
                products.append(s)  # keep original if parsing fails

        '''
        # Pick main product (largest heavy atom count)
        def main_product(smiles_list):
            return max(smiles_list, key=lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms())

        main_prod = main_product(products)

        # Add steps for each reactant
        for react in reactants:
            all_steps.append((react, ", ".join(reagents), main_prod))
        '''
        
        def product_size(smiles):
            mol = Chem.MolFromSmiles(smiles)
            return mol.GetNumHeavyAtoms() if mol else -1

        main_prod = max(products, key=product_size)

        # Add one step per reactant
        for react in reactants:
            mol = Chem.MolFromSmiles(react)
            if mol:
                react = Chem.MolToSmiles(mol, canonical=True)
            all_steps.append((react, reagents, main_prod))
            # all_steps.append((react, ", ".join(reagents), main_prod))

    return all_steps

# -------------------------------
# Drawing function
# -------------------------------
def draw_multistep_pathway(steps, output_path):
    """
   #Draws a clean multistep pathway with molecules, arrows, and reagents.
   #`steps`: list of tuples (reactant, reagent, product)
    """
    if not steps:
        raise ValueError("No steps to draw.")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Extract stage molecules and reagents
    stage_smiles = [steps[0][0]] + [p for (_, _, p) in steps]
    reagents = [r for (_, r, _) in steps]

    mols = [Chem.MolFromSmiles(s) for s in stage_smiles]
    for mol in mols:
        Chem.rdDepictor.Compute2DCoords(mol)

    subImgSize = (250, 250)
    mol_imgs = [Draw.MolToImage(mol, size=subImgSize) for mol in mols]

    arrow_width = 120
    canvas_width = len(mols) * subImgSize[0] + (len(mols)-1) * arrow_width
    canvas_height = subImgSize[1] + 80

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    try:
        font = ImageFont.truetype("arial.ttf", 20)
    except:
        font = ImageFont.load_default()

    x = 0
    for i, mol_img in enumerate(mol_imgs):
        canvas.paste(mol_img, (x, 50))
        x += subImgSize[0]

        if i < len(reagents):
            arrow_start = x
            arrow_end = x + arrow_width - 20
            y = subImgSize[1] // 2 + 50

            draw.line((arrow_start, y, arrow_end, y), fill="black", width=3)
            draw.polygon([(arrow_end, y), (arrow_end-10, y-5), (arrow_end-10, y+5)], fill="black")

            bbox = draw.textbbox((0, 0), reagents[i], font=font)
            text_w = bbox[2] - bbox[0]
            text_h = bbox[3] - bbox[1]
            text_x = (arrow_start + arrow_end) // 2 - text_w // 2
            text_y = y - text_h - 10
            draw.text((text_x, text_y), reagents[i], font=font, fill="black")

            x += arrow_width

    canvas.save(output_path)
    print(f"Pathway saved as {output_path}")

# -------------------------------
# Example usage
# -------------------------------
if __name__ == "__main__":
    smirks_example = "CCO>H2SO4>C=C;C=C>Br2>CCBr;CCBr>NH3>CCN"
    # smirks_example = "C1=CC=CC=C1.O=C(O)CCl>[Cu]Br2.CN(C)C.CS(=O)(=O)[O-].[K+]>C1=CC=CC=C1C(=O)NC"
    steps = parse_multistep_smirks(smirks_example)
    output_file = os.path.join("data/reactions", "Multi_step_pathway.png")
    draw_multistep_pathway(steps, output_file)


