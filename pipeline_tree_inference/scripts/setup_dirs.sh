#!/bin/bash
#SBATCH -c 1                    # number of "cores"
#SBATCH -t 0-01:00:00           # time in d-hh:mm:ss
#SBATCH -p general              # partition
#SBATCH -q public               # QOS
#SBATCH -e slurm.%A_%a.err      # file to save STDERR for each sub-job
#SBATCH --export=NONE           # Purge the job-submitting shell environment
#SBATCH --mail-type=FAIL
#SBATCH --mail-user="%u@asu.edu"

# Create all directories necessary for analyses
echo "Creating directories"

# 1. Group/lab directory --for output of all pipeline steps
# As a member of Qiyun"s lab, you should have access to the lab directory

# Parameters
ks=(10 50 100 200 400 600 800 1000)
ps=(0)
directories=('alignments' 'trim' 'trees' 'select_criteria/presence_absence_copies')
database="${1:?ERROR -- must pass a database name}"
path=/data/qzhu44/marker_genes_project/${USER}/wol2/${database}/copy_number_first_strategy/min_marker_genes_per_genome/bit_score_threshold_1.0

for k in ${ks[@]};do
    for p in ${ps[@]}; do
        for dir in ${directories[@]};
            do mkdir -p ${path}/k_${k}_p_${p}/${dir}
        done
    done
done


echo "Done"


