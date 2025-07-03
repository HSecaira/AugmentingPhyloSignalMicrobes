#!/usr/bin/env python3
'''Script to parse the output of HMMer
'''

# Imports
import os, argparse, re
import numpy as np
import pandas as pd

################################################################################################
# Functions

def parse_hmmer_outputs(genomes, path):
    results = []
    for genome in genomes:
        # Load file
        with open(f'{path}/{genome}_hits.txt', 'r') as f:
            for line in f:
                if not line.startswith('#') and len(line.strip()) != 0:
                    tmp = re.sub('\s+', ' ', line.strip()).split(' ')
                    # Ignore elements of description columns
                    tmp = tmp[:18]
                    tmp.append(genome)
                    results.append(tmp)

    # Convert to dataframe
    df = pd.DataFrame(results, columns = ['target', 'acession_1', 'query', 'acession_2', 'e-value_full', 'score_full', 'bias_full',
                        'e-value_domain', 'score_domain', 'bias_domain',
                        'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'genome'])
    df.set_index('target', inplace = True)
    df = df.astype({'score_full': 'float64'})
    return df

def get_orfs(df):
    orfs = {}
    genes = df['query'].unique()
    for scg in genes:
        orfs[scg] = []
        filtered_df = df[df['query'] == scg]
        grouped = filtered_df.groupby('genome')
        for genome, group in grouped:
            if not group.empty:
                # Take only the best hit, according to bit score
                top_idx = group.sort_values(by = 'score_full', ascending = False).index[0]
                orfs[scg].append(top_idx)
    return orfs

def save_orfs(orfs, path):
    for scg, ofs in orfs.items():
        with open(f'{path}/{scg}.txt', 'w') as f:
            for of in ofs:
                f.write(f'{of}\n')

def load_proteins(genomes, path):
    proteins = {}
    for genome in genomes:
        with open(f'{path}/{genome}.faa', 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    proteins[header] = ''
                else:
                    proteins[header] += line.strip()

    return proteins

def save_proteins(orfs, proteins, path):
    for scg, ofs in orfs.items():
        with open(f'{path}/{scg}.fasta', 'w') as f:
            for of in ofs:
                f.write(f'>{of}\n{proteins[of]}\n')

################################################################################################

# # Parser
parser = argparse.ArgumentParser()
parser.add_argument('--markers', type = str, required = True)
args = parser.parse_args()

################################################################################################

# Parameters
replicates = np.arange(0, 10, dtype = int)
# Paths
genomes_path = f''

for rep in replicates:
    print(f'At replicate: {rep}')
    # Load genomes
    genomes_path = f'/mnt/scratch0/hsecaira/BipNet/hmm/data/wol2/taxonomic_imbalance/out/replicate_{rep}/'
    genomes = np.loadtxt(f'{genomes_path}/genomes/genomes.txt', dtype = str)
    # Parse HMMer outputs
    dataPathOut = f'/mnt/scratch0/hsecaira/BipNet/hmm/{args.markers}/taxonomic_imbalance/replicate_{rep}'
    df = parse_hmmer_outputs(genomes, f'{dataPathOut}/hits')
    # Get ORFs
    orfs = get_orfs(df)
    print(f'Number of markers: {len(orfs)}')
    # Save ORFs
    save_orfs(orfs, f'{dataPathOut}/orfs')
    # Load proteins
    proteins = load_proteins(genomes, f'{genomes_path}/input_genomes')
    # Save proteins of ORFs
    save_proteins(orfs, proteins, f'{dataPathOut}/proteins')
    # Save markers
    markers = df['query'].unique()
    np.savetxt(f'{dataPathOut}/all_genes.txt', markers, fmt = '%s')
