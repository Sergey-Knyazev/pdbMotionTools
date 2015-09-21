#!/usr/bin/env bash

TEMP=`getopt --long -o "o:f:d:h" "$@"`
eval set -- "$TEMP"
SCRIPT_BASE=`basename $0`
SCRIPT_DIR=`dirname $0`

SPLIT_STATES_PATH="$SCRIPT_DIR/../pdbSidechainsBuilder/pdbSplitStates.sh"

TEMP_DIR=".temp"
PDB_MODELS_DIR="pdb_models"

function print_help () {
    echo "$SCRIPT_BASE - calculate free energy of protein conformations"
    echo " "
    echo "$SCRIPT_BASE [options]"
    echo " "
    echo "options:"
    echo "-h        show brief help"
    echo "-f        specify the pdb file"
    echo "-d        specify the directory with pdb files"
    echo "-o        specify out csv file"
}

while true ; do
    case "$1" in
        -h )
            print_help
            exit 0
        ;;
        -f )
            IN_FILE=$(readlink -f $2)
            shift 2
        ;;
        -d )
            IN_DIR=$(readlink -f $2)
            shift 2
        ;;
        -o )
            OUT_CSV=$(readlink -f $2)
            shift 2
        ;;
        *)
            break
        ;;
    esac
done

if ( [ -n "$IN_FILE" ]  && [ -n "$IN_DIR" ] ) || ( [ -z "$IN_FILE" ]  && [ -z "$IN_DIR" ] ); then
    echo "Specify either pdb file or directory with pdb files but not both!"
    print_help
    exit 1
fi

if [ -z "$OUT_CSV" ]; then
    echo "Out csv file should be specified!"
    print_help
    exit 1
fi

if [ -n "$IN_FILE" ] && [ ! -f "$IN_FILE" ]; then
    echo "Input pdb file is not found!"
    print_help
    exit 1
fi

if [ -n "$IN_DIR" ] && [ ! -d "$IN_DIR" ]; then
    echo "Directory with pdb files is not found!"
    print_help
    exit 1
fi

if [ -e "$TEMP_DIR" ]; then
    i="0"
    while true ; do
        TEMP_DIR=".temp"$i
        if [ ! -e "$TEMP_DIR" ]; then
            break
        fi
        i=$[$i+1]
    done
fi

mkdir "$TEMP_DIR"
cd "$TEMP_DIR"
mkdir "$PDB_MODELS_DIR"
if [ -n "$IN_FILE" ]; then
    "$SPLIT_STATES_PATH" "$IN_FILE" "$PDB_MODELS_DIR"
    IN_DIR=$TEMP_DIR
fi

for f in $PDB_MODELS_DIR/*.pdb
do
	filename=$(basename "$f")
	filename="${filename%.*}"
    mkdir "$filename"
    cd "$filename"
	log_file=${filename}.gromacs.log

	echo $(basename "$f") processing........

	echo --------------------------------------- >> $log_file
	echo $(basename "$f") processing >> $log_file

    pdb2gmx -f "../$f" -p ${filename}.top -o ${filename}.gro -ff gromos53a6 -water none -ignh &>> $log_file
	grompp -f $SCRIPT_DIR/energy_check.mdp -c ${filename}.gro -p ${filename}.top \
		 -o ${filename}.tpr &>> $log_file
    editconf -f ${filename}.gro -o ${filename}_gro.pdb &>> $log_file
    mdrun -s ${filename}.tpr -rerun ${filename}_gro.pdb &>> $log_file
    cd ..
done

for f in $PDB_MODELS_DIR/*.pdb
do
	filename=$(basename "$f")
	filename="${filename%.*}"
    t=$(echo $(sed -n '/Kinetic/{n;p}' $filename/md.log) | cut -d " " -f1)
    echo "$filename;$t" >> $OUT_CSV
done

rm -rf "$TEMP_DIR"
cd "$SCRIPT_DIR"
