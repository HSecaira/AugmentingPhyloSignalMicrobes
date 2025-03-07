conda activate eggnog-mapper
exedir=/home/drz/Programs/EggNOG-mapper/2.1.7
export PATH=$exedir:$exedir/eggnogmapper/bin:$PATH
cp -r /mnt/store0/dbs/eggnog /dev/shm/
export EGGNOG_DATA_DIR=/dev/shm/eggnog
tmpdir=$(mktemp -d)
query=/mnt/store0/qiyun/wol2/data/concat.faa
emapper.py -i $query -o output --cpu 32 --temp_dir $tmpdir
rm -rf $tmpdir
rm -rf /dev/shm/eggnog
