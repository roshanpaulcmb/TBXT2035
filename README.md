# TBXT2035
https://tbxtchallenge.org/

# Overview
This project's goal is to find druggable binding pockets for TBXT, an intrinsically disordered protein.

1. Generating predicted structures using AlphaFold3
2. Running (enhanced sampling) MD simulations
3. Finding druggable binding pockets
4. Docking pharmocophores to binding pokects
5. Running (enhanced sampling) MD simulations with lead ligand candidates
6. Reiteration

# MD Simulations conducted

* Implicit solvent
* Explicit solvent
* REST2 explicit solvent

# Tools

* dcdTools - for editing trajectory files
* analysisTools - simulation analysis and visualization

# Background

TBXT (Brachyury) is a T-box transcription factor that is essential for notochord formation and early embryonic patterning.
TBXT has a structured DNA- binding region and a disordered C-terminal IDR.
TBXT binds to DNA as a dimer and its G177D mutant is implicated in chordoma in humans.
Currently, TBXT is deemed undruggable.
We hypothesize that molecular dynamics (MD) simulations can reveal transient cryptic pockets suitable for drug screening.
