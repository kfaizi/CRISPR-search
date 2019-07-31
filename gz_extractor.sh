#!/bin/bash
# gz_extractor.sh
# extract NCBI WGS fasta files of form *.fsa_nt.gz

src=$1 # directory with gzipped fastas
dest=$2 # directory for unzipped fastas

mkdir -vp $dest

cd $src

for i in $( ls *.fsa_nt.gz ); do # for every zipped fasta in src,
    zipped=$i
    gunzip -kc < $zipped > $dest/${i%%.gz} # unzip a copy and send to dest file of same name; keeps src 
done

exit
