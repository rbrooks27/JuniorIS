from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image, ImageDraw, ImageFont


# Creating a molecule from SMILES
m = Chem.MolFromSmiles("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2")
img = Draw.MolToImage(m, size=(600, 600))  # You can adjust size if you want
png_file = "data/testing.png"
img.save(png_file)
print(f"Molecule image saved as '{png_file}'")

# Creating a reaction from SMARTS
rxn_smarts = "CC(=O)O.CCO>[H+]>CC(=O)OCC.O"
rxn = AllChem.ReactionFromSmarts(rxn_smarts)
img = Draw.ReactionToImage(rxn)
png_file = "data/testing_reaction.png"
img.save(png_file)
print(f"Reaction image saved as '{png_file}'")

# Creating a 3 step synthesis scheme image
rxn_smarts = "CC(=O)O.CCO>[H+]>CC(=O)OCC.O"
rxn = AllChem.ReactionFromSmarts(rxn_smarts)
rxn_smarts2 = "CC(=O)OCC.O.CCO>[H+]>CC(=O)OCCOCC.O"
rxn2 = AllChem.ReactionFromSmarts(rxn_smarts2)
rxn_smarts3 = "CC(=O)OCCOCC.O.CCO>[H+]>CC(=O)OCCOCCOCC.O"
rxn3 = AllChem.ReactionFromSmarts(rxn_smarts3)
img1 = Draw.ReactionToImage(rxn)
img2 = Draw.ReactionToImage(rxn2)
img3 = Draw.ReactionToImage(rxn3)
# Combine images horizontally
total_width = img1.width + img2.width + img3.width
max_height = max(img1.height, img2.height, img3.height)
combined_img = Image.new('RGB', (total_width, max_height), (255, 255, 255))
combined_img.paste(img1, (0, 0))
combined_img.paste(img2, (img1.width, 0))
combined_img.paste(img3, (img1.width + img2.width, 0))
png_file = "data/testing_synthesis_scheme.png"
combined_img.save(png_file)
print(f"Synthesis scheme image saved as '{png_file}'")

# Creating a multi-reaction scheme image with annotations
rxn_smarts_a = "CC(=O)O.CCO>[H+]>CC(=O)OCC.O"
rxn_a = AllChem.ReactionFromSmarts(rxn_smarts_a)
rxn_smarts_b = "CC(=O)OCC.O.CCO>[H+]>CC(=O)OCCOCC.O"
rxn_b = AllChem.ReactionFromSmarts(rxn_smarts_b)
img_a = Draw.ReactionToImage(rxn_a)
img_b = Draw.ReactionToImage(rxn_b)
# Combine images vertically with annotations
spacing = 50  # Space for annotation
total_height = img_a.height + img_b.height + spacing
max_width = max(img_a.width, img_b.width)
combined_img2 = Image.new('RGB', (max_width, total_height), (255, 255, 255))
combined_img2.paste(img_a, (0, 0))
combined_img2.paste(img_b, (0, img_a.height + spacing))
# Add annotations
draw = ImageDraw.Draw(combined_img2)
font = ImageFont.load_default()
draw.text((10, img_a.height + 10), "Step 1: Esterification", fill="black", font=font)
draw.text((10, img_a.height + spacing + 10), "Step 2: Chain Extension", fill="black", font=font)
png_file = "data/testing_multi_reaction_scheme.png"
combined_img2.save(png_file)
print(f"Multi-reaction scheme image saved as '{png_file}'")


# Creating a multi-step reaction scheme image with arrows, reagents,  and reactants/products names underneath the molecules.
rxn_smarts1 = "CC(=O)O.CCO>[H+]>CC(=O)OCC.O"
rxn1 = AllChem.ReactionFromSmarts(rxn_smarts1)
rxn_smarts2 = "CC(=O)OCC.O.CCO>[H+]>CC(=O)OCCOCC.O"
rxn2 = AllChem.ReactionFromSmarts(rxn_smarts2)
img1 = Draw.ReactionToImage(rxn1)
img2 = Draw.ReactionToImage(rxn2)
# Combine images horizontally with arrows and annotations
arrow = Image.new('RGB', (50, max(img1.height, img2.height)),   (255, 255, 255))
draw = ImageDraw.Draw(arrow)
draw.polygon([(10, arrow.height // 2 - 10), (40, arrow.height // 2), (10, arrow.height // 2 + 10)], fill="black")
total_width = img1.width + arrow.width + img2.width
max_height = max(img1.height, img2.height)
combined_img3 = Image.new('RGB', (total_width, max_height + 50), (255, 255, 255))
combined_img3.paste(img1, (0, 0))
combined_img3.paste(arrow, (img1.width, 0))
combined_img3.paste(img2, (img1.width + arrow.width, 0))
# Add annotations
draw = ImageDraw.Draw(combined_img3)
png_file = "data/testing_multi_step_reaction_scheme.png"
combined_img3.save(png_file)
print(f"Multi-step reaction scheme image saved as '{png_file}'")




