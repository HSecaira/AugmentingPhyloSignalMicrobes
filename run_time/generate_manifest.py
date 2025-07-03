#!/usr/bin/env python3
'''Script to generate manifest file for running
jobs in parallel
'''

# Imports
import os
import numpy as np

################################################################################################

# Parameters
replicates = np.arange(0, 10, dtype = int)
Ss = np.logspace(1, 4.7, num = 7, dtype = int)
Gs = np.logspace(1.7, 5, num = 7, dtype = int)
# Variables
path_script = f'/mnt/store0/BipNet/new_analysis/scripts/run_time'

with open('manifest.txt', 'w') as f:
    # Iterate over parameters
    for rep in replicates:
        for G in Gs:
            for S in Ss:
                # Number of markers to select
                if G <= 10000:
                    ks = np.logspace(1, round(np.log10(G), 2), num = 7, dtype = int)
                    # For the last value, substract one
                    ks[-1] = ks[-1] - 1
                else:
                    ks = np.logspace(1, 4, num = 7, dtype = int) # up to 10,000
                for k in ks:
                    # Generate command
                    cmd = f'/usr/bin/time -v python -u {path_script}/select_markers.py --S {S} --G {G} --k {k} --rep {rep}'
                    f.write(f'{cmd}\n')
