#!/usr/bin/env python3
'''Script to generate species and gene trees
from simulated data of a presence-absence matrix
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

def greedy_power_mean_sample_final(data, k, p, pseudocount, min_universality_genes):
    """Select k rows from a matrix such that the selection criterion by column is maximized.

    Parameters
    ----------
    data : ndarray (2D)
        Input data matrix.
    k : int
        Number of rows to select.
    p : float
        Exponent.

    Returns
    -------
    ndarray (1D)
        Row indices selected in order.
    """

    n, m = data.shape

    # Matrix is empty
    if n == 0 or m == 0:
        raise ValueError(f'Matrix is empty!')

    # Matrix contains only zeroes
    if (data == 0).all():
        raise ValueError(f'Matrix only contains 0\'s')

    if (data.sum(axis = 1) == min_universality_genes).all():
        raise ValueError(f'Matrix only contains genes present in less than {min_universality_genes}')

    if k >= n:
        raise ValueError(f'k should be smaller than {n}')
    
    # Add pseudocount
    data = data + pseudocount

    # cumulative gene counts
    counts = np.zeros(m, dtype = int)

    # gene indices in original data matrix
    indices = np.arange(n)

    # indices of selected genes
    selected = []

    # will select k genes sequentially
    for i in range(k):
        # calculate counts after adding each gene
        sums_ = counts + data

        # select a gene that maximizes the power mean gene count per genome, using the cumulative matric
        if isinstance(p, int) or isinstance(p, np.int64): 
            choice = pmean(sums_, int(p), axis = 1).argmax()
        elif p == 'min':
            choice = sums_.min(axis = 1).argmax()
        elif p == 'max':
            choice = sums_.max(axis = 1).argmax()
        else:
            raise ValueError(f'Invalid p: {p}.')

        # append index of selected gene
        selected.append(indices[choice])

        # update per-species gene counts
        counts = sums_[choice]

        # remove selected gene from data matrix
        data = np.delete(data, choice, axis = 0)

        # remove selected gene from indices
        indices = np.delete(indices, choice)

    return np.sort(np.array(selected))

def add_noise_trees(trees, percentage_tips = 0.5, seed = 42):
    '''Function that add tree arrangements, specifically the tips are shuffled
    using a random generator.
    Parameters
    ----------
    trees : List of skbio TreeNodes
        Gene trees
    percentage_tips : float
        Percentage of tips to shuffle
    Returns
    -------
    List of skbio TreeNodes
        Gene trees with tips shuffled
    '''
    shuffled_trees = []
    # rev = lambda items: items.reverse()
    rng = np.random.default_rng(seed = seed)

    for tree in trees:
        tree_tmp = tree.copy()
        n_tips = tree_tmp.count(tips = True)
        k = int(n_tips * percentage_tips)
        if k > 2:
            shuffled_trees.append([st for st in tree_tmp.shuffle(k = k, shuffle_f = rng.shuffle)][0])
        else:
            shuffled_trees.append(tree_tmp)

    return shuffled_trees

def marker_gene_selection(data):
    '''Function that selects marker genes based on k (number of genes to select)
        and p (exponent of power mean --> selection criterion).
    Parameters
    ----------
    data : Numpy 2D matrix of ints
        Gene/rows species/columns matrix with entries between 0 and max_copies
    Returns
    -------
    Dictionary
        Keys : combination of k and p
        Values : index of selected genes
    '''

    # Parameters
    n, m = data.shape
    # ks = np.linspace(10, n - 1, dtype = int, num = 9) # Number of marker genes
    ks = [100]
    ps = [-100, -75, -50, -25, 0, 25, 50, 75, 100] # Exponent of power mean

    # Marker gene selection
    marker_genes = {}
    for i, k in enumerate(ks):
        for j, p in enumerate(ps):
            # Select marker genes
            select = greedy_power_mean_sample_final(data, k = k, p = p, pseudocount = 0.1, min_universality_genes = 1)
            # Add to results
            marker_genes[f'select_k_{k}_p_{p}'] = select

    return marker_genes

def save_marker_gene_trees(gene_trees, rep, shuffled, path, noise, k, p):
    '''Function that saves the gene trees of selected marker genes into a file.
    Parameters
    ----------
    gene_trees : List of skbio TreeNodes
        Gene trees
    marker_genes :     Dictionary
        Keys : combination of k and p
        Values : index of selected genes
    rep : int
        Replicate number
    shuffled : Boolean
        Indicates whether trees are shuffled
    path : str
        Output path
    Returns
    -------
    None
    '''
    if shuffled:
        name = f'{path}/marker_genes_shuffled'
    else:
        name = f'{path}/marker_genes'

    with open(f'{name}_select_k_{k}_p_{p}_noise_{noise}_rep_{rep}.nwk', 'w') as f:
        for tree in gene_trees:
            # Rename nodes
            for node in tree.tips():
                node.name = node.name.split(' ')[0]
            f.write(f'{str(tree)}')

################################################################################################

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('--noise', type = float, required = True)
parser.add_argument('--p', type = int, required = True)
parser.add_argument('--k', type = int, required = True)
parser.add_argument('--rep', type = int, required = True)
args = parser.parse_args()

################################################################################################

# replicates = np.arange(1, 51, dtype = int)
dataPathIn = '/mnt/store0/BipNet/new_analysis/tools/zombi-code/ZOMBI/Tutorials/Tutorial2/500_replicates'

print(f'Replicate: {args.rep}')
# Load species tree
tree = TreeNode.read(f'{dataPathIn}/replicate_{args.rep}/output/T/ExtantTree.nwk')
# Load edges
edges_families, edges_genomes = load_gene_families(tree, f'{dataPathIn}/replicate_{args.rep}/output/G/Genomes')
# Build presence-absence matrix
adj, genomes, ogs = adjMatv2(edges_families, edges_genomes)
# Reformat
adj1 = adj
adj2 = adj1.copy()
adj1[adj1 > 1] = 1
# Copy presence/absence matrix
adj_mod = adj1.copy()
# Select genes to remove
remove = remove_genes(adj_mod)
# Remove genes
adj_mod = np.delete(adj_mod, remove, axis = 0)
# Remove genes from list
ogs_mod = np.delete(ogs, remove, axis = 0)
# Copy matrix
adj_mod2 = adj2.copy()
# Remove genes
adj_mod2 = np.delete(adj_mod2, remove, axis = 0)
print(f'\tPresence-absence matrix shape: {adj_mod.shape}')
print(f'\tParameters: noise = {args.noise}; p = {args.p}; k = {args.k}')
# Select marker genes
select = greedy_power_mean_sample_final(adj_mod, k = args.k, p = args.p, pseudocount = 0.1, min_universality_genes = 1)
marker_genes = [ogs_mod[i] for i in select]
# Get gene trees
gene_trees = [TreeNode.read(f'{dataPathIn}/replicate_{args.rep}/output/G/Gene_trees/{og}_prunedtree.nwk') for og in marker_genes]
# Add noise to gene trees
shuffled_trees = add_noise_trees(gene_trees, percentage_tips = args.noise, seed = args.rep)
# Save marker genes
path = f'{dataPathIn}/replicate_{args.rep}/marker_genes'
save_marker_gene_trees(gene_trees, rep = args.rep, shuffled = False, path = path, noise = args.noise, k = args.k, p = args.p)
save_marker_gene_trees(shuffled_trees, rep = args.rep, shuffled = True, path = path, noise = args.noise, k = args.k, p = args.p)

print(f'Done!')



