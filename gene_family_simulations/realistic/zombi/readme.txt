This directory contains the scripts necessary to simulate
a species tree and gene families along its branches using
Zombi simulator (https://academic.oup.com/bioinformatics/article/36/4/1286/5578480).

Each simulation was repeated for 25 replicates, fixing a seed for replicability.
Parameters for species trees anf genome were taken from Tutorial 2 (see Wiki: https://github.com/AADavin/Zombi/wiki); see also methods section of our paper.

After generating commands for simulation (generate_files_cmmands.py), run them in parallel with:
parallel -j 2 < replicate_x/commands/commands.txt &> run_x.log & disown
# Remove temperary files of Zombi simulations otherwise, they can fill up the entire
# space in your machine
find . -type d -name "Gene_families" -exec rm -r {} +
find . -type f -name "*completetree*" -exec rm {} +

---
We aslo provide script used to select markers for each replicate.
