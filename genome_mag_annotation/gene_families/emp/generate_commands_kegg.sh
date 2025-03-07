#!/bin/bash

# Data paths
input=/mnt/store1/data/emp500/mags/prokka/faa
output=/mnt/scratch0/hsecaira/BipNet/emp500/kegg/raw

# Read bin
while read bin; do
	echo "bash ${output}/kofamscan.run ${input}/${bin}.faa ${output}/out/${bin}.tsv ${output}/profiles/prokaryote.hal ${output}/ko_list" >> ${output}/commands_kegg.txt
done < ${output}/bins.txt
