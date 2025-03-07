#!/usr/bin/env python3
'''Script to generate species and gene trees
from simulated data of a copy number matrix.
Usage: python special_cases_multi_copy.py --max_copies 1 --noise 0.25
'''

# Imports
import os, argparse
import numpy as np
from scipy.stats import pmean
from skbio import TreeNode
from skbio.tree import nj
from skbio.stats.distance import DistanceMatrix

################################################################################################
# Functions

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

def jaccard_distance(si, sj):
    '''Jaccard similarity index.
    Parameters
    -----------
    si, sj : sets
            Containing genes of species i and j (str or int)
    Returns
    -------
    jaccard_index : float
                    Between 0 and 1.
    '''

    intersection = len(si.intersection(sj))
    union = len(si.union(sj))
    if union > 0:
        dist = 1 - (intersection / union)
    else:
        dist = 1.0

    return dist

def get_gene_trees(tree, data, sps):
    '''Function that returns all the gene trees with the same topology as the species tree
        but some species are sheared according to the pattern of the presence absence matrix.
    Parameters
    ----------
    tree : scikit-bio TreeNode
        Species tree
    data : numpy 2D array; int
        Presence absence matrix. Rows indicate genes and columns species
    sps : numpy 1D array; str
        Species names
    Returns
    -------
    List of scikit-bio TreeNodes
        Gene trees
    '''
    gene_trees = []
    for i in range(len(data)):
        tree_copy = tree.copy()
        num_sps = len(data[i, :].nonzero()[0])
        # Avoid empty genes and ORFans
        if num_sps > 1:
            gt = tree_copy.shear(sps[data[i, :].nonzero()])
            gt.prune()
        else:
            gt = TreeNode()
        gene_trees.append(gt)
    return gene_trees

def add_copies_tree(trees, data, sps):
    '''Function that returns a gene tree with the same topology as the species tree
        with species removed based on the presence absence matrix. Then, gene copies
        are added as sister branches based on the copy numbers.
    Parameters
    ----------
    trees : List of scikit-bio TreeNodes
        Gene trees
    data : numpy 2D array; int
        Presence absence matrix with copy numbers. Rows indicate genes and columns species
    sps : numpy 1D array; str
        Species names
    Returns
    -------
    List of Scikit-bio TreeNodes
        Gene trees
    '''
    multicopy_gene_trees = []
    # Iterate over genes
    for i, tree in enumerate(trees):
        tree_ = tree.copy()
        gene = data[i]
        # Iterate over gene in species idx
        for idx in gene.nonzero()[0]:
            num_copies = gene[idx]
            # Avoid single-copy genes and ORFans
            if num_copies > 1 and len(gene.nonzero()[0]) > 1:
                sp = sps[idx]
                node = tree_.find(sp)
                # Iterate over gene copies
                for copy in range(num_copies - 1):
                    leaf = TreeNode(name = f'{sp}_{copy}', length = node.length + node.median)
                    # Add copy as branch sister
                    node.parent.append(leaf)
        # Add tree to list
        multicopy_gene_trees.append(tree_)

    return multicopy_gene_trees

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

def nodes_median_depths(trees):
    '''Function that calculates the median and depth of nodes
        of trees.
    Parameters
    ----------
    trees : List of skbio TreeNodes
        Gene trees
    Returns
    -------
    List of skbio TreeNodes
        Gene trees with median and depths of nodes
    '''
    for tree in trees:
        if tree.count(tips = True) > 0:
            for node in tree.postorder(include_self = True):
                if node.length is None:
                    node.length = 0.0
                if node.is_tip():
                    node.depths = [0.0]
                    node.median = 0.0
                else:
                    node.depths = [y + x.length for x in node.children for y in x.depths]
                    node.median = np.median(node.depths)
    return trees

def simulate_gene_species_matrix(G, S, seed = 42, max_copies = 1):
    '''Function that simulates a presence absence matrix of G genes
        over S species and then adds copy numbers.
    Parameters
    ----------
    G : int
        Number of gene families
    S : int
        Number of species
    seed : int
        Seed for reproducibility
    max_copies : int
        Maximum number of copies of a gene family in a species
    Returns
    -------
    Numpy 2D matrix of ints
        Gene/rows species/columns matrix with entries between 0 and max_copies
    '''
    rng = np.random.default_rng(seed = seed)

    # Presence-absence matrix
    fraction_genes = rng.uniform(low = 0.1, high = 1.0, size = S)
    data = [rng.choice([0, 1], size = G, p = [1 - x, x]) for x in fraction_genes]
    data = np.vstack(data).T
    data_ = data.copy()
    n, m = data_.shape
    labels = np.array([f'{label}' for label in range(m)])
    genes = np.arange(n, dtype = int)

    # Gene copy numbers over matrix
    for i in range(n):
        gene = data_[i, :]
        # Sample copies
        copies = rng.integers(low = 1, high = max_copies, size = len(gene.nonzero()[0]), endpoint = True)
        # Add to data_
        gene[gene > 0] = copies

    return data_, labels, genes

def build_sps_tree(labels, genes, data):
    '''Function that builds a species tree based on the
        gene-species matrix with Jaccard distance and 
        Neighbor Joining.
    Parameters
    ----------
    labels : numpy 1D array of ints
        dsadas
    genes : numpy 1D array of ints or str
        Gene names
    data : Numpy 2D matrix of ints
        Gene/rows species/columns matrix with entries between 0 and max_copies
    Returns
    -------
    skbio TreeNode
        Species tree
    '''
    n, m = data.shape
    # Calculate distance between species
    dm = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            if j > i:
                dm[i, j] = dm[j, i] = jaccard_distance(set(genes[data[:, i].nonzero()]), set(genes[data[:, j].nonzero()]))

    # Neighbor joining
    dm_skbio = DistanceMatrix(dm, labels)
    sps_tree_nj = nj(dm_skbio)
    
    return sps_tree_nj

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
    ks = np.linspace(1, n - 1, dtype = int, num = 11) # Number of marker genes
    ps = [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100] # Exponent of power mean

    # Marker gene selection
    marker_genes = {}
    for i, k in enumerate(ks):
        for j, p in enumerate(ps):
            # Select marker genes
            select = greedy_power_mean_sample_final(data, k = k, p = p, pseudocount = 0.1, min_universality_genes = 1)
            # Add to results
            marker_genes[f'select_k_{k}_p_{p}'] = select

    return marker_genes

def save_marker_gene_trees(gene_trees, marker_genes, rep, shuffled, path):
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
        name = f'{path}/marker_genes/marker_genes_shuffled'
    else:
        name = f'{path}/marker_genes/marker_genes'

    for key, v in marker_genes.items():
        with open(f'{name}_{key}_rep_{rep}.nwk', 'w') as f:
            for gene in v:
                tree = gene_trees[gene]
                # Rename nodes
                for node in tree.tips():
                    node.name = node.name.split('_')[0]
                f.write(f'{str(tree)}')

################################################################################################

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('--max_copies', type = int, required = True)
parser.add_argument('--noise', type = float, required = True)
args = parser.parse_args()

################################################################################################

replicates = np.arange(0, 25, dtype = int)
G = 50
S = 10
for rep in replicates:
    print(f'Replicate: {rep}')
    # Simulate data. G and S must be somewhat larger.
    data, labels, genes = simulate_gene_species_matrix(G = G, S = S, seed = rep, max_copies = args.max_copies)
    # Build species tree
    sps_tree_nj = build_sps_tree(labels, genes, data)
    # Save trees
    path = f'./simulations_G_{G}_S_{S}/multicopy_trees_{args.max_copies}_copy_noise_{args.noise}'
    os.system(f'mkdir -p {path}/marker_genes')
    os.system(f'mkdir -p {path}/sps_trees')
    os.system(f'mkdir -p {path}/pipeline')
    sps_tree_nj.write(f'{path}/sps_trees/sps_tree_rep_{rep}.nwk')
    # Gene trees
    gene_trees = get_gene_trees(sps_tree_nj, data, labels)
    # Calculate median and node depths
    gene_trees = nodes_median_depths(gene_trees)
    # Add gene copies as sister branches
    multicopy_gene_trees = add_copies_tree(gene_trees, data, labels)
    # Add noise to gene trees.
    shuffled_multicopy_gene_trees = add_noise_trees(multicopy_gene_trees, percentage_tips = args.noise, seed = rep)
    # Select marker genes
    marker_genes = marker_gene_selection(data)
    # Save marker gene trees
    save_marker_gene_trees(multicopy_gene_trees, marker_genes, rep, shuffled = False, path = path)
    save_marker_gene_trees(shuffled_multicopy_gene_trees, marker_genes, rep, shuffled = True, path = path)


print(f'Done!')



