 ## 3D Protein Folding Simulation
## Overview

This project simulates the process of protein folding in a 3D space using computational techniques. It uses a simplified HP (Hydrophobic-Polar) model and applies Monte Carlo simulation with Simulated Annealing to find a low-energy protein structure.

## Objective

The main goal is to arrange amino acids in a way that minimizes the total energy of the system, especially by bringing hydrophobic (H) residues closer together, similar to real biological protein folding.

## Features
3D lattice-based protein model
Monte Carlo random movement simulation
Simulated Annealing optimization
Energy minimization technique
3D visualization of protein structure
Energy vs Steps graph
Temperature vs Steps graph
Auto-save output image
## Technologies Used
Python
NumPy
Matplotlib
## How to Run
Install required libraries:
pip install numpy matplotlib
Run the script:
python protein_folding_3d.py
Output:
3D visualization window
Saved image in outputs/protein_folding_result.png
## Model Used

The project uses the HP model, where:

H (Hydrophobic) → tends to cluster together
P (Polar) → interacts with surroundings

Energy is calculated based on non-bonded H-H contacts.

## Results

The simulation starts with a random structure and gradually improves it using probabilistic optimization. The final output shows a folded structure with reduced energy along with graphs for better understanding.

## Future Improvements
Real protein data (PDB files) integration
Advanced moves (crankshaft, pivot)
Interactive 3D visualization
Animation of folding process
Better optimization techniques
## Author

Sunil Matwa
25BCE10466
