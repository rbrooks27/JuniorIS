
"""
reactionGraph.py

Builds a reaction graph from (reactant, reagents, product) tuples.

Features:
- RDKit canonicalization to avoid mismatches
- Forward adjacency list  : reactant → [(product, reagents)]
- Reverse adjacency list  : product → [(reactant, reagents)]
- Breadth First Search (BFS) retrosynthesis
- Forward exploration for planning
- Optional deduplication
- Optional exporting later (JSON, gpickle, etc.)

This file is intended to integrate with combined123.py
which will take a list of steps and draw a multistep synthesis PNG.
"""

'''
from collections import defaultdict, deque # This helps with graph structures specifically for adjacency lists and BFS/DFS


class ReactionGraph:
    def __init__(self):
        # adjacency list: reactant -> list of (product, reagents)
        self.graph = defaultdict(list)
        
        # reverse graph for retrosynthesis: product -> reactants
        self.reverse_graph = defaultdict(list)
    def add_reaction(self, reactant, reagents, product):
        # Add a reaction step to the directed graph
        self.graph[reactant].append((product, reagents))
        self.reverse_graph[product].append((reactant, reagents))
    def load_steps(self, steps):
        # Steps: list of (reactant, reagents, product) tuples
        for react, reagents, prod in steps:
            self.add_reaction(react, reagents, prod)
    def find_path(self, start_materials, target):
        # Use BFS to find the shortes synthesis route from start_materials to target
        start_materials = set(start_materials) # for quick lookup
        queue = deque([(target, [])]) # (current molecule, path)
        visited = set([target]) # to prevent cycles
        while queue:
            molecule, path = queue.popleft() # current molecule and path to it
            # Check if we've reached a starting material
            if molecule in start_materials: 
                return list(reversed(path)) # found a path, return reversed
            
            # Look backwards: find reactions that produce this molecule
            for reactant, reagents in self.reverse_graph.get(molecule, []):
                # This conditional check if we've already visited this reactant and if not add to queue
                if reactant not in visited:
                    visited.add(reactant)
                    new_step = (reactant, reagents, molecule)
                    queue.append((reactant, path + [new_step]))
        return None # no path found
'''

from collections import defaultdict, deque
from rdkit import Chem


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
def canonical(smiles: str) -> str:
    """Convert SMILES to canonical RDKit form."""
    if not smiles:
        return smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles  # preserve invalid strings
    return Chem.MolToSmiles(mol, canonical=True)



# ---------------------------------------------------------
# Reaction Graph
# ---------------------------------------------------------
class ReactionGraph:
    def __init__(self):
        self.graph = defaultdict(list)          # reactant → [(product, reagents)]
        self.reverse_graph = defaultdict(list)  # product → [(reactant, reagents)]

    # -----------------------------------------------------
    # Add a individual reaction
    # -----------------------------------------------------
    def add_reaction(self, r, reag, p):
        r = canonical(r)
        p = canonical(p)

        self.graph[r].append((p, reag))
        self.reverse_graph[p].append((r, reag))

    # -----------------------------------------------------
    # Load dataset including semicolon-multistep lines
    # -----------------------------------------------------
    def load_dataset(self, path):
        steps = []

        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # split multi-step lines: r>reag>p; r>reag>p; ...
                chunks = line.split(";")

                for chunk in chunks:
                    chunk = chunk.strip()
                    if not chunk:
                        continue

                    parts = chunk.split(">")
                    if len(parts) != 3:
                        print(f"Skipping malformed: {chunk}")
                        continue

                    r, reag, p = parts
                    r, reag, p = r.strip(), reag.strip(), p.strip()

                    steps.append((r, reag, p))

        # actually store them in graphs
        self.load_steps(steps)
        return steps

    def load_steps(self, steps):
        for r, reag, p in steps:
            self.add_reaction(r, reag, p)

    # -----------------------------------------------------
    # Retrosynthesis BFS search
    # -----------------------------------------------------
    def find_path(self, start_materials, target):
        start_materials = {canonical(x) for x in start_materials}
        target = canonical(target)
        

        queue = deque([(target, [])])
        visited = set([target])

        while queue:
            molecule, path = queue.popleft()

            # success condition
            if molecule in start_materials:
                return list(reversed(path))

            # backward search
            for reactant, reagents in self.reverse_graph.get(molecule, []):
                if reactant not in visited:
                    visited.add(reactant)
                    new_step = (reactant, reagents, molecule)
                    queue.append((reactant, path + [new_step]))

        return None

    # -----------------------------------------------------
    # Export entire pathway as single SMIRKS string
    # -----------------------------------------------------
    def path_to_smirks(self, path):
        if not path:
            return None
        smirks_steps = [f"{r}>{reag}>{p}" for (r, reag, p) in path]
        return ";".join(smirks_steps)


# -----------------------------
# Basic Test Case
# -----------------------------
if __name__ == "__main__":
    """
    print("\n=== ReactionGraph Dataset Test ===\n")

    G = ReactionGraph()
    G.load_dataset("dataset/textbookReactions.txt")

    start = ["CCO"]
    target = "COCC"

    print("Finding retrosynthesis path...\n")
    path = G.find_path(start, target)

    if not path:
        print("No pathway found.")
    else:
        for r, reag, p in path:
            print(f"{r} --[{reag}]→ {p}")

        print("\nSMIRKS:")
        print(G.path_to_smirks(path))
    """

    """
    # Example: 3-step benzylamine → acetamide derivative synthesis
    test_steps = [
        ("c1ccccc1", "CuBr2", "BrCCc1ccccc1",
         "c1ccccc1>CuBr2>BrCCc1ccccc1"),
        ("BrCCc1ccccc1", "NH3", "NCc1ccccc1",
         "BrCCc1ccccc1>NH3>NCc1ccccc1"),
        ("NCc1ccccc1", "AcCl", "O=CNCCc1ccccc1",
         "NCc1ccccc1>AcCl>O=CNCCc1ccccc1")
    ]

    start_materials = ["c1ccccc1"]
    target = "O=CNCCc1ccccc1"

    G = ReactionGraph()
    G.load_steps(test_steps)

    print("\nFinding path...")
    path = G.find_path(start_materials, target)

    if path is None:
        print("No path found.")
    else:
        print("\n=== Retrosynthetic Path Found ===")
        for react, reag, prod, smirks in path:
            print(f"{react} --[{reag}]→ {prod}")

        combined = G.path_to_smirks(path)

        print("\n=== Combined SMIRKS ===")
        print(combined)

        print("\n=== Test Complete ===")
"""




