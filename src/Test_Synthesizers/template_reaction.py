"""
Template for visualizing chemical reactions using RDKit.
- Define the reaction using SMARTS/SMIRKS.
- Generate 2D coordinates for reactants and products.
- Kekulize molecules to make bond types explicit.
- Save the reaction diagram as a PNG.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os

# -------------------------------
# Setup
# -------------------------------
os.makedirs("data", exist_ok=True)

# -------------------------------
# 1. Define the reaction SMARTS
# -------------------------------
rxn_smarts = "C1(=CC=C(C=C1C[P](=O)(OCC)OCC)[S](F)(F)(F)(F)F)[N+](=O)[O-].C1(=CC=CC=C1)C=O>N(C)(C)C=O.CC(C)(C)[O-].[K+]>C1(=CC=C(C=C1C=CC2=CC=CC=C2)[S](F)(F)(F)(F)F)[N+](=O)[O-]"
rxn = AllChem.ReactionFromSmarts(rxn_smarts)

rxn2_smarts = "CC(=O)O.CCO>[H+]>CC(=O)OCC.O"
rxn2 = AllChem.ReactionFromSmarts(rxn2_smarts)

# -------------------------------
# 2. Generate 2D coordinates for reactants/products
# -------------------------------
reactants = list(rxn.GetReactants())
products = list(rxn.GetProducts())

reactants2 = list(rxn2.GetReactants())
products2 = list(rxn2.GetProducts())

for mol in reactants + products:
    if mol is None:
        continue
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass
    AllChem.Compute2DCoords(mol)

for mol in reactants2 + products2:
    if mol is None:
        continue
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass
    AllChem.Compute2DCoords(mol)
# -------------------------------
# 3. Draw and save the reaction
# -------------------------------
# ReactionToImage returns a PIL.Image object
img = Draw.ReactionToImage(rxn)
png_file = "data/reaction1.png"
img.save(png_file)
print(f"Reaction image saved as '{png_file}'")

img = Draw.ReactionToImage(rxn2)
png_file = "data/reaction2.png"
img.save(png_file)
print(f"Reaction image saved as '{png_file}'")

# -------------------------------
# 4. Optional: save individual reactants/products
# -------------------------------
'''
for idx, mol in enumerate(reactants, start=1):
    Draw.MolToFile(mol, f"data/reactant_{idx}.png")
for idx, mol in enumerate(products, start=1):
    Draw.MolToFile(mol, f"data/product_{idx}.png")
print("Individual reactant/product images saved in 'data/'")
'''
