#!/usr/bin/env python3
'''Script to simulate copy number matrices
'''

# Imports
import os
import numpy as np
import argparse
from scipy.stats import pmean

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
    # labels = np.array([f'{label}' for label in range(m)])
    # genes = np.arange(n, dtype = int)

    # Gene copy numbers over matrix
    for i in range(n):
        gene = data_[i, :]
        # Sample copies
        copies = rng.integers(low = 1, high = max_copies, size = len(gene.nonzero()[0]), endpoint = True)
        # Add to data_
        gene[gene > 0] = copies

    return data_

################################################################################################

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('--max_copies', type = int, required = True)
parser.add_argument('--outdir', type = str, required = True)
args = parser.parse_args()

################################################################################################

replicates = np.arange(0, 10, dtype = int)
Ss = np.logspace(1, 4.7, num = 7, dtype = int)
Gs = np.logspace(1.7, 5, num = 7, dtype = int)
# Gs = [50, 100, 200]
# Ss = [10, 50, 100]
print(f'Simulating data for {len(replicates)} replicates, {len(Gs)} gene families and {len(Ss)} species')
for replicate in replicates:
    for G in Gs:
        for S in Ss:
            print(f'\tAt replicate: {replicate}, G: {G}, S: {S}')
            # Simulate copy number matrix
            data = simulate_gene_species_matrix(G, S, seed = replicate, max_copies = args.max_copies)
            # Save matrix as binary file with memmap
            filename = f'{args.outdir}/matrix_S_{S}_G_{G}_rep_{replicate}.dat'
            fp = np.memmap(filename, dtype = 'int32', mode = 'w+', shape = data.shape)
            fp[:] = data[:]
            fp.flush()

print(f'Done!')

