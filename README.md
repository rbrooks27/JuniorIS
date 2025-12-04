# OrgoNizer - (Multistep Synthesis Organic Chemistry)
The purpose of this project is to design and implement an interactive software system that teaches multistep synthesis through decision trees and constraint satisfaction techniques. By breaking synthesis into a sequence of guided yes/no decisions, the system will reinforce fundamentals and provide educational feedback at each step. 
While professional chemistry tools exist—such as ChemDraw or Reaxys--these platforms are not designed for novice learners. This project fills the gap by prioritizing transparency, interpretability, and pedagogy over raw prediction power.

## Features
    - Loads a curated set of textbook organic reactions  
    - Canonicalizes molecule and reaction SMILES for consistency  
    - Builds a directed graph (**nodes = molecules**, **edges = reactions**)  
    - Generates **single-step and multistep reaction pathways** using BFS  
    - Presents synthesis steps in a transparent, educational format  
    - Supports exploration of alternative synthetic routes  
    - Designed specifically for **learning**, not prediction

## Requirements & Setup

### Update your Python
```pip install python==3.10.18```

### Install, create, activate a Conda environment
```conda install -c conda-forge rdkit=2022.09.5 pillow```
```conda create -n rdkit-env python=3.10.18```
```conda activate rdkit-env```

### How to exit Conda environment
```conda deactivate```

## How to Run

### Navigate to the src/ directory
```cd src```

### Configure inputs inside generatePathways.py

#### Modification of Dataset File
Change the ```dataset_file``` variable

#### Max Number of Dataset Lines
Change the ```MAX_LINES``` variable

### Target + Starting Materials
    ``` bash
    target_raw = " "   # target molecule (SMILES)
    starting_materials_raw = ["", ""]  # list of starting materials
    ```

### Run the Pathway Generator
```python generatePathways.py```


## Repository Structure
```bash
├── README.md
├── main.pdf                 # Full academic write-up (Junior IS)
├── src/                     # Source code for the OrgoNizer tool
│   ├── reactionGraph.py     # Reaction graph builder + search engine
│   ├── generatePathways.py  # Pathway generation script (BFS)
│   ├── combined123.py       # Utility and testing functions
│   ├── data/
│   │   ├── reactions/       # Reaction Generated Images
│   │
│   ├── datasets/
│   │       ├── reactionSmilesFigShare2024.txt   # Reaction SMILES 1.37 million 
│   │       ├── textbookReactions.txt            # Personally curated textbook reactions
│   
└── ...
```


## Feature Calendar

| **Issue** | **Due date** | |
| --------- | ------------ | -- ||
| [Example issue description with link](https://github.com/hmm34/example-annotated-bibliography/issues/1) | 9/1/25 |Completed|
| [Find Dataset for Organic Molecule based on Attributes needed](https://github.com/rbrooks27/JuniorIS/issues/1) | 9/19/25 |Completed|
| [Building database to store all molecules and reactions (if needed)](https://github.com/hmm34/example-annotated-bibliography/issues/1) | 9/19/25 |Was Not Needed|
| [Develop Decision Tree Logic](https://github.com/rbrooks27/JuniorIS/issues/3) | 11/1/25 |Completed|
| [Input molecule by name or drawing](https://github.com/rbrooks27/JuniorIS/issues/4) | 0/0/00 |Completed|
| [Confirm and display input molecules](https://github.com/rbrooks27/JuniorIS/issues/5) | 0/0/00 |Completed|
| [Search function for molecules](https://github.com/rbrooks27/JuniorIS/issues/6) | 0/0/00 |Completed|
| [Generate single-step transformations](https://github.com/rbrooks27/JuniorIS/issues/7) | 0/0/00 |Completed|
| [Constraint enforcement on pathways.](https://github.com/rbrooks27/JuniorIS/issues/8) | 11/1/25|Completed|
| [Generate multistep synthesis pathways](https://github.com/rbrooks27/JuniorIS/issues/9) | 0/0/00 |Completed|
| [User interface for interaction](https://github.com/rbrooks27/JuniorIS/issues/10) | 0/0/00 |Completed|
| [Step-by-step walkthrough mode with explanations](https://github.com/rbrooks27/JuniorIS/issues/11) | 0/0/00 |Completed|
| [Test on actual synthesis problems](https://github.com/rbrooks27/JuniorIS/issues/12) | 11/1/25|Completed|