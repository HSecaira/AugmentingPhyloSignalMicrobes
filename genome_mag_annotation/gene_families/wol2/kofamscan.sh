query=/mnt/store0/qiyun/wol2/data/concat.faa
cwd=$(pwd)
output=$cwd/output.tsv
cd /home/drz/Programs/KofamScan/1.3.0
$(which time) ./exec_annotation -f detail-tsv --no-report-unannotated -o $output $query
cd $cwd
