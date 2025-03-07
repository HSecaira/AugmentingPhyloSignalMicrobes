#!/usr/bin/env python3

'''
Extract mapping between contig and prokka ids
'''
# Imports
import os

#################################################################################################

def parse_gff(file):
    '''
    Parse GFF file and extract the contig and prokka id
    '''
    contigs = {}
    lengths = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('##sequence-region'):
                l = line.strip().split(' ')
                lengths[l[1]] = int(l[-1])
            else:
                l = line.strip().split('\t')
                if len(l) == 9 and l[8].startswith('ID='):
                    contig_id = l[0]
                    prokka_id = l[8].split(';')[0].split('=')[-1]

                    if contig_id not in contigs:    
                        contigs[contig_id] = [prokka_id]
                    else:
                        contigs[contig_id].append(prokka_id)
    return contigs, lengths


def save_mapping(file, contigs):
	'''
	Save mapping of contigs and prokka ids into a file
	'''
	with open(file, 'w') as f:
		for contig, idxs in contigs.items():
			for idx in idxs:
				f.write(f'{contig},{idx}\n')

def save_lengths(file, lengths):
	with open(file, 'w') as f:
		for contig, length in lengths.items():
			f.write(f'{contig},{length}\n')

def load_file(file):
	data = []
	with open(file, 'r') as f:
		for line in f:
			data.append(line.strip())
	return data

#################################################################################################

# Paths
pathIn = '/mnt/store1/data/emp500/mags/prokka/gff'

# Load bin names
bins = load_file('./bins.txt')

# Iterate
for bi in bins:
	print(f'Bin: {bi}')
	# Parse gff
	contigs, lengths = parse_gff(f'{pathIn}/{bi}.gff')
	# Save mapping
	save_mapping(f'./{bi}.mapping', contigs)
	# Save lengths
	save_lengths(f'./{bi}.lengths', lengths)

print(f'Done!')
