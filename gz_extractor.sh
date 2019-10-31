#!/bin/bash
# gz_extractor.sh

src=$1 # directory with gzipped fastas
dest=$2 # directory for unzipped fastas

mkdir -vp $dest

cd $src

for i in $( ls *.gz ); do # for every zipped file in src,
    zipped=$i
    gunzip -kc < $zipped > $dest/${i%%.gz} # unzip to dest; keeps src copy
done

exit
