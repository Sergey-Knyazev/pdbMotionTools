#!/usr/bin/env bash

TEMP=`getopt --long -o "i:o:h" "$@"`
eval set -- "$TEMP"
SCRIPT_BASE=`basename $0`
SCRIPT_DIR=`dirname $0`
while true ; do
    case "$1" in
        -h )
                echo "$SCRIPT_BASE - optimizes transformation"
                echo " "
                echo "$SCRIPT_BASE [options] application [arguments]"
                echo " "
                echo "options:"
                echo "-h        show brief help"
                echo "-i        specify a directory with transformation"
                echo "-o        specify a directory to store output in"
                exit 0
        ;;
        -i )
            IN_DIR=$2
            shift 2
        ;;
        -o )
            OUT_DIR=$2
            shift 2
        ;;
        *)
            break
        ;;
    esac
done

for f in $IN_DIR/*.pdb
do
	filename=$(basename "$f")
	filename="${filename%.*}"
	pdb2gmx -f $f -p $OUT_DIR/${filename}.top -o $OUT_DIR/${filename}.gro -ff gromos53a6 -water none -ignh
	grompp -f $SCRIPT_DIR/min.mdp -c $OUT_DIR/${filename}.gro -p $OUT_DIR/${filename}.top -o $OUT_DIR/${filename}_input_min.tpr
	mdrun -s $OUT_DIR/${filename}_input_min.tpr -deffnm $OUT_DIR/${filename}_min -v
	editconf -f $OUT_DIR/${filename}_min.gro -o $OUT_DIR/${filename}_min.pdb
done


