#!/usr/bin/env bash

# minimazes the energy of pdb structures by GROMACS

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
                echo "-i        specify a directory wich contain pdb files for energy minimization"
                echo "-o        specify an output directory"
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

if [ "$(ls -A $OUT_DIR)" ]; then
    echo Error: output folder is not empty!
    exit 1
fi

for f in $IN_DIR/*.pdb
do
	filename=$(basename "$f")
	filename="${filename%.*}"
	log_file=$OUT_DIR/${filename::-4}.gromacs.log
	energy_log=$OUT_DIR/${filename::-4}.energy_relax.log

	echo $(basename "$f") processing........

	echo --------------------------------------- >> $log_file
	echo $(basename "$f") processing >> $log_file

	pdb2gmx -f $f -p $OUT_DIR/${filename}.top -o $OUT_DIR/${filename}.gro -ff gromos53a6 -water none -ignh &>> $log_file
	grompp -f $SCRIPT_DIR/min.mdp -c $OUT_DIR/${filename}.gro -p $OUT_DIR/${filename}.top \
		 -o $OUT_DIR/${filename}_input_min.tpr &>> $log_file
	mdrun -s $OUT_DIR/${filename}_input_min.tpr -deffnm $OUT_DIR/${filename}_min -v &>> $log_file
	editconf -f $OUT_DIR/${filename}_min.gro -o $OUT_DIR/${filename}_min.pdb &>> $log_file
done

cat $log_file | grep -e Step -e "processing$" >> $energy_log

