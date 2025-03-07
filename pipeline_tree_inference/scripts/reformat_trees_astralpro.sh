#!/bin/bash
#SBATCH -c 1                    # number of "cores"
#SBATCH -t 0-01:00:00           # time in d-hh:mm:ss
#SBATCH -p general              # partition
#SBATCH -q grp_lliu80		# QOS
#SBATCH -e slurm.%A_%a.err      # file to save STDERR for each sub-job
#SBATCH --export=NONE           # Purge the job-submitting shell environment
#SBATCH --mail-type=FAIL
#SBATCH --mail-user="%u@asu.edu"

# Manifest file is provided
# Execute script as: sbatch reformat_trees_astralpro.sh manifest

# Set up sbatch parameters and environment
module purge

# Reformat with sed
ps=(0)
ks=(10 50 100 200 400 600 800 1000)

# Data paths
# Remember to create directory and subdirectories
th_copy=1.0
database="${1:?ERROR -- must pass a database name}"
input=/data/qzhu44/marker_genes_project/${USER}/wol2/${database}/copy_number_first_strategy/min_marker_genes_per_genome/bit_score_threshold_${th_copy}

for k in ${ks[@]}; do
    for p in ${ps[@]}; do
        echo "Criterion: k_${k}_p_${p}"
        sed -E 's/([A-Z]{8})_[a-z0-9]{1,5}/\1/g' ${input}/k_${k}_p_${p}/select_criteria/presence_absence_copies/select_k_${k}_p_${p}_shrink.trees > ${input}/k_${k}_p_${p}/select_criteria/presence_absence_copies/select_k_${k}_p_${p}_multitree.nwk
    done
done

