
'''
reactionGraph.py
- reading cleaned reaction steps
- canonicalizing molecules
- building adjacency list {reactant: [(product, reagents)]}
- preventing duplicates and cycles
- exporting structure (json / .gpickle / text / svg later)

The goal here is to use Directed Acyclic Graph (DAG) or Djikstra-like structures 
to represent reaction pathways for efficient searching and visualization. I will use 
molecules/immediate steps as nodes and reactions/reagents as edges.
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