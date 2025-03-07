Annotation of ORFs
---
* EggNog

eggnog-mapper v.2.1.12 installed via conda
        `mamba install -c bioconda -c conda-forge eggnog-mapper`
on environment eggnog-mapper.

eggNOG release 5.0 was used as the reference for annotation.
        Path: /mnt/store0/dbs/eggnog

# Run 4 jobs in parallel
parallel -j 4 < commands_eggnog.txt &> run.log & disown

---
* KEGG

kofamscan v1.3.0 installed via conda
        `mamba install bioconda::kofamscan`
on environment kofamscan

Kofam release 2022-03-01 was used as the reference for annotation.
        Path: /mnt/store0/dbs/kegg/kofam/2022-03-01     

# Run 4 jobs in parallel
parallel -j 4 < commands_eggnog.txt &> run.log & disown

# In output from Kofamscan, the best hit for every KO is higlighted with an asterisk
# The best hit refers to a hit with a score higher than the threshold  defined for the KO.
# See https://github.com/takaram/kofam_scan for details.
# Command to select only the best hits
sed -n 's/^\*[^\t]*\t//p' file > output 

