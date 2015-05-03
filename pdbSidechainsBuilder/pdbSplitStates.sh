#!/usr/bin/env bash

# Separates the multi-modal pdb file to many single-model pdb files.
# usage: ./pdbSplit states input.pdb out/patn/

filename=$(basename "$1")
filename="${filename%.*}"
command="%dp $1 > $2/${filename}_%03d.pdb\n"

grep -n 'MODEL\|ENDMDL' $1 | cut -d: -f 1 | \
 awk -v c="$command" '{if(NR%2) printf "sed -n %d,",$1+1; else printf c , $1-1,NR/2;}' |  bash -sf
