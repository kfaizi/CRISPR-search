#!/bin/bash
# gz_extractor.sh
# helps methodically extract fasta from NCBI WGS gz files,
# of form "*.fsa_nt.gz"

# destination directory must exist??


src=$1 # pass the name of directory w/ .gzs
dest=$2 # pass the name of directory for fastas

mkdir -v -p $dest

cd $src

for i in $( ls *.fsa_nt.gz ); do # for every zipped fasta in src,
    zipped=$i
    gunzip -kc < $zipped > $dest/${i%%.gz} # unzip a copy and send to dest file of same name; keeps src 
done

exit
