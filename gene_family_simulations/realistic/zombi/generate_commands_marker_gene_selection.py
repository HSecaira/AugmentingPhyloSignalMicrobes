#!/usr/bin/env python3
'''Script to generate commands for
scriipt `marker_gene_selection_600.py`
'''

# Imports
import os
import numpy as np
import argparse
from scipy.stats import pmean
from skbio import TreeNode

################################################################################################
# Functions

def load_gene_families(tree, path):
    '''
    Function that loads the gene families of ZOMBI simulations
    as an edge list
    Parameters
    ----------
        tree: skbio TreeNode containing species
        path: str containing the path to the directory that has genomic information
    Returns
    -------
        edges: Arrays containing the genome and gene family
    '''
    genomes, families = [], []
    for node in tree.tips():
        file = os.path.join(path, f'{node.name}_GENOME.tsv')
        with open(file, 'r') as f:
            for line in f:
                if line[0] != 'P':
                    gene_family = line.strip().split('\t')[1]
                    genomes.append(node.name)
                    families.append(gene_family)
    return np.array(families), np.array(genomes)

def adjMatv2(edges_ogs, edges_genomes):
    # Get unique elements
    ogs, ogs_indices = np.unique(edges_ogs, return_inverse = True)
    genomes, genomes_indices = np.unique(edges_genomes, return_inverse = True)
    # Calculate bin counts for each combination of indices
    counts = np.bincount(ogs_indices * len(genomes) + genomes_indices, minlength = len(ogs) * len(genomes))
    # Reshape counts as adjacency matrix
    adj = counts.reshape(len(ogs), len(genomes))
    return adj, genomes, ogs

def remove_genes(adj_mod):
    '''
    Function that removes genes if:
        present in less than 4 genomes/species
    '''
    remove = np.array([i for i in range(len(adj_mod)) if adj_mod[i].sum() < 4])
    return remove

################################################################################################

# Parser
# parser = argparse.ArgumentParser()
# parser.add_argument('--noise', type = float, required = True)
# args = parser.parse_args()

################################################################################################

dataPathIn = '/mnt/store0/BipNet/new_analysis/tools/zombi-code/ZOMBI/Tutorials/Tutorial2/500_replicates'

#replicates = np.arange(24, 51, dtype = int)
replicates = [1, 5, 6, 8, 11]
for rep in replicates:
    # Load species tree
    tree = TreeNode.read(f'{dataPathIn}/replicate_{rep}/output/T/ExtantTree.nwk')
    # Load edges
    edges_families, edges_genomes = load_gene_families(tree, f'{dataPathIn}/replicate_{rep}/output/G/Genomes')
    # Build presence-absence matrix
    adj, genomes, ogs = adjMatv2(edges_families, edges_genomes)
    # Reformat
    adj1 = adj
    adj1[adj1 > 1] = 1
    # Copy presence/absence matrix
    adj_mod = adj1.copy()
    # Select genes to remove
    remove = remove_genes(adj_mod)
    # Remove genes
    adj_mod = np.delete(adj_mod, remove, axis = 0)
    # Remove genes from list
    ogs_mod = np.delete(ogs, remove, axis = 0)
    # Parameters for marker gene selection
    n, m = adj_mod.shape
    ks = np.geomspace(10, n - 1, dtype = int, num = 9) # Number of marker genes
    ps = [-100, -75, -50, -25, 0, 25, 50, 75, 100] # Exponent of power mean
    noises = [0.05, 0.1, 0.25, 0.5]
    # Write commands
    os.system(f'mkdir -p {dataPathIn}/replicate_{rep}/commands')
    os.system(f'mkdir -p {dataPathIn}/replicate_{rep}/marker_genes')
    with open(f'{dataPathIn}/replicate_{rep}/commands/commands.txt', 'w') as f:
        for noise in noises:
            for k in ks:
                for p in ps:
                    f.write(f'python3 -u {dataPathIn}/marker_gene_selection.py --noise {noise} --p {p} --k {k} --rep {rep}\n')

print(f'Done!')



