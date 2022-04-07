#!/bin/bash
# gz_extractor.sh

src=$1 # directory with zipped fastas
dest=$2 # directory for unzipped fastas

mkdir -vp $dest

cd $src

# for gzips
for i in $( ls *.gz ); do # for every zipped file in src,
    gzipped=$i
    gunzip -kc < $gzipped > $dest/${i%%.gz} # unzip to dest; keeps src copy
done

# for zips
for i in $( ls *.zip ); do
    zipped=$i
    unzip -d $zipped $dest/${i%%.zip}
done

exit
