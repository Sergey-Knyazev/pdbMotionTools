#!/usr/bin/env bash

OXYGEN_ADDER_DIR=oxygenAdder
STATE_SPLITTER_DIR=''

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
                echo "-i        specify a file with transformation"
                echo "-o        specify a directory to store output in"
                exit 0
        ;;
        -i )
            IN_FILE=$2
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

xbase=${IN_FILE##*/}

$SCRIPT_DIR/$OXYGEN_ADDER_DIR/oxygen_adder -i $IN_FILE -o $OUT_DIR/$xbase
$SCRIPT_DIR/$STATE_SPLITTER_DIR/pdbSplitStates.sh $OUT_DIR/$xbase $OUT_DIR/
rm $OUT_DIR/$xbase
for f in $OUT_DIR/*.pdb
do
#    Scwrl4 -i $f -o $f >> ${f%.*}.log
    Scwrl4 -i $f -o $f
done
