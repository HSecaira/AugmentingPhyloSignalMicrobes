#!/usr/bin/env python3
'''Script to select k markers from a copy number matrix
'''

# Imports
import os
import numpy as np
import argparse
from scipy.stats import pmean

################################################################################################
# Functions

def greedy_power_mean_sample_final(data, k, p, pseudocount):
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

    return np.array(selected)

def greedy_power_mean_sample_mask(data, k, p, pseudocount = 0.1):
    """Iteratively select rows from a matrix such that the selection criterion by
    column is maximized.

    Parameters
    ----------
    data : ndarray of shape (2D)
        Input data matrix (gene by species).
    k : int
        Number of rows (genes) to select.
    p : float or str
        Exponent of generalized mean, or special values: "min" or "max".
    pseudocount : float, optional
        Pseudocount to add to each cell value.
    min_universality_genes : int, optional
        Minimum universality threshold.

    Returns
    -------
    ndarray (1D)
        Row indices selected in order.

    """
    
    # numbers of genes (n) and species (m)
    n, m = data.shape

    # matrix is empty
    if n == 0 or m == 0:
        raise ValueError('Matrix is empty!')

    # matrix contains only zeroes
    if (data == 0).all():
        raise ValueError('Matrix only contains 0\'s')

    if k >= n:
        raise ValueError(f'k should be smaller than {n}')

    # Add pseudocount
    data = data.astype(np.float64, copy = False) + pseudocount

    # cumulative gene counts
    counts = np.zeros(m, dtype = np.float64)

    # gene indices in original data matrix
    all_original_indices = np.arange(n)
    # Boolean mask to keep track of available genes
    available_row_mask = np.ones(n, dtype = bool)

    # indices of selected genes
    selected = np.empty(k, dtype = int)

    # iteratively select k genes
    for i in range(k):

        # Create new array containing only available rows
        current_data_slice = data[available_row_mask, :]
        # calculate counts after adding each gene
        sums_ = counts + current_data_slice

        # select a gene that maximizes the power mean gene count per species, using the cumulative matrix
        if isinstance(p, int) or isinstance(p, np.int64): 
            choice = pmean(sums_, int(p), axis = 1).argmax()
        elif p == 'min':
            choice = sums_.min(axis = 1).argmax()
        elif p == 'max':
            choice = sums_.max(axis = 1).argmax()
        else:
            raise ValueError(f'Invalid p: {p}.')

        # append index of selected gene
        original_indices_of_available_rows = all_original_indices[available_row_mask]
        chosen_original_idx = original_indices_of_available_rows[choice]
        selected[i] = chosen_original_idx

        # update per-species gene counts
        counts = sums_[choice, :]

        # Update mask
        available_row_mask[chosen_original_idx] = False

    # return np.sort(np.array(selected))
    return selected

################################################################################################

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('--G', type = int, required = True)
parser.add_argument('--S', type = int, required = True)
parser.add_argument('--rep', type = int, required = True)
parser.add_argument('--k', type = int, required = True)
args = parser.parse_args()

################################################################################################

# Load copy number matrix
print(f'Selecting {args.k} genes from a matrix with {args.G} genes and {args.S} species, rep {args.rep}...')
data = data = np.memmap(f'/mnt/store0/BipNet/new_analysis/scripts/run_time/data/matrix_S_{args.S}_G_{args.G}_rep_{args.rep}.dat',
                dtype = 'int32', mode = 'r', shape = (args.G, args.S))

# Select k genes as markers
# selected_genes = greedy_power_mean_sample_final(data = data, k = args.k, p = 0, pseudocount = 0.1)
selected_genes = greedy_power_mean_sample_mask(data = data, k = args.k, p = 0, pseudocount = 0.1)
print(f'Done!, selected genes: {len(selected_genes)}')

