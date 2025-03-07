#!/usr/bin/env python3
'''Script to generate SpeciesTreeParameters.tsv and GenomeParameter.tsv files
and commands for Zombi simulations using realistic speciation, extinction, 
duplication, transfer, and loss rates.
'''


import os
import numpy as np

####################################################################################

# Software Path
zombi_dir = f'/mnt/store0/BipNet/new_analysis/tools/zombi-code/ZOMBI'

# Path variables
path = f'/mnt/store0/BipNet/new_analysis/tools/zombi-code/ZOMBI/Tutorials/Tutorial2/500_replicates'

# Parameters for simulation
replicates = np.arange(1, 25, dtype = int)

with open(f'{path}/commands_species_trees.txt', 'w') as ft, open(f'{path}/commands_genomes.txt', 'w') as fg:
    for rep in replicates:
        print(f'Replicate {rep}')
        # Create directories
        os.system(f'mkdir -p {path}/replicate_{rep}/output')
        # Replace seed in parameter files
        os.system(f'sed -e "s/SEED 42/SEED {rep}/" {path}/SpeciesTreeParameters.tsv > {path}/replicate_{rep}/SpeciesTreeParameters_rep_{rep}.tsv')
        os.system(f'sed -e "s/SEED 42/SEED {rep}/" {path}/GenomeParameters.tsv > {path}/replicate_{rep}/GenomeParameters_rep_{rep}.tsv')
        # Create command files
        ft.write(f'python3 {zombi_dir}/Zombi.py T {path}/replicate_{rep}/SpeciesTreeParameters_rep_{rep}.tsv {path}/replicate_{rep}/output\n')
        fg.write(f'python3 {zombi_dir}/Zombi.py Gm {path}/replicate_{rep}/GenomeParameters_rep_{rep}.tsv {path}/replicate_{rep}/output\n')


print(f'Done!')




