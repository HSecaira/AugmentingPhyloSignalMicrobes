#!/bin/bash

# Data paths
input=/mnt/store1/data/emp500/mags/prokka/faa
output=/mnt/scratch0/hsecaira/BipNet/emp500/eggnog/raw

# Read bin
while read bin; do
	echo "bash ${output}/eggnog-mapper.run ${input}/${bin}.faa ${output}/out/${bin}" >> ${output}/commands_eggnog.txt
done < ${output}/bins.txt
