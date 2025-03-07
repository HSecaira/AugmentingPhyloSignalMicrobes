# Running phylophlan v3.1.68

Remember to activate the mamba environment
	mamba activate phylophlan

# By default, phylophlan discards genomes with less than 100 or 34 marker genes if the database if phylophlan or amphora, respectively.

1. Input genomes
	Input genomic data is located in ./emp500|wol2/input_genomes with a soft link to the original path.

2. Generate configuration file
	phylophlan_write_config_file \
	-d a \
	-o data.cfg \
	--db_aa diamond \
	--map_aa diamond \
	--msa mafft \
	--trim trimal \
	--tree1 iqtree \
	--verbose 2>&1 | tee phylophlan_write_config_file.log

3. Reconstruct the tree of life. Change -d depending on the databse (phylophlan or amphora).
	phylophlan \
	-i input_genomes \
	-d phylophlan \
	-f data.cfg \
	--diversity high \
	--fast \
	-o output \
	--nproc 4 \
	--verbose 2>&1 | tee phylophlan.log

# Phylophlan generates a tmp directory which contains the marker genes (what we want).

# Changing the minimum number of marker genes per genome

1. Input genomes
        Input genomic data is located in ./emp500|wol2/input_genomes with a soft link to the original path.

2. Generate configuration file
        phylophlan_write_config_file \
        -d a \
        -o data.cfg \
        --db_aa diamond \
        --map_aa diamond \
        --msa mafft \
        --trim trimal \
        --tree1 iqtree \
        --verbose 2>&1 | tee phylophlan_write_config_file.log

3. Reconstruct the tree of life. Change -d depending on the databse (phylophlan or amphora).
        phylophlan \
        -i input_genomes \
        -d phylophlan \
        -f data.cfg \
        --diversity high \
	--min_num_markers 1 \
        --fast \
        -o output \
        --nproc 10 \
        --verbose 2>&1 | tee phylophlan.log


