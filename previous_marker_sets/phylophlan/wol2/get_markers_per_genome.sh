# Go to phylophlan output directory
# output/tmp/markers_aa
for f in $(ls .); do echo -e "$(echo $f | cut -d "." -f 1):$(less $f | grep -c ">")"; done > markers_per_genome.txt

