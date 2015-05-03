#!/usr/bin/env bash

# scripts implements following procedures:
# - separates multi-model pdb file to many single-model pdb;
# - adds sidechains to protein backbone for every single-model pdb by using Sqwrl4 tool.

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
                echo "-i        specify an input pdb file"
                echo "-o        specify an output directory"
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

if [ "$(ls -A $OUT_DIR)" ]; then
    echo "Error: the output folder is not empty!"
    exit 1
fi

xbase=${IN_FILE##*/}
log_file=$OUT_DIR/${xbase%.*}.Scwrl4.log

$SCRIPT_DIR/$OXYGEN_ADDER_DIR/oxygen_adder -i $IN_FILE -o $OUT_DIR/$xbase
$SCRIPT_DIR/$STATE_SPLITTER_DIR/pdbSplitStates.sh $OUT_DIR/$xbase $OUT_DIR/
rm $OUT_DIR/$xbase
for f in $OUT_DIR/*.pdb
do
    echo $(basename "$f") pricessing....
    echo --------------------------------------- >> $log_file
    echo $(basename "$f") processing >> $log_file
    Scwrl4 -i $f -o $f >> $log_file
done
