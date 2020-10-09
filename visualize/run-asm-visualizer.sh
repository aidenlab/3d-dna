#!/bin/bash

#### Description: Script to visualize draft assemblies: a wrapper around remap-contig-mnd-to-asm-mnd.awk and Juicebox pre. Also dumps 2D annotation files.
#### Usage: bash ./run-asm-visualizer.sh [options] path_to_cprops path_to_asm <contig-mnd-file>.
#### Input: cprops file, asm file, mnd file.
#### Output: hic file, scaffold and super-scaffold 2D annotation files. Optional: contig annotation file?
#### Parameters: zoom (default is >=1 calculated to fit the assembly chromosome); mapq threshold (0, 1, 30, all; default is 1 for mapq>=1).
#### Options: -l <path_to_gap_bed>, -p <true/false> to use parallelization for speed-up.
#### Dependencies: Java; GNU Parallel if available.
#### Written by: Olga Dudchenko & Sanjit Batra, version date 07/18/2016

USAGE="
*****************************************************
Visualizing draft genomes in juicebox: 18 July 2016

USAGE: ./run-asm-visualizer.sh [options] <path_to_cprops> <path_to_asm> <path_to_contig_mnd> 

DESCRIPTION:
This is a script to visualize draft assemblies (represented in input by their cprops and asm files) from pairwise contact data represented by Juicer merged_nodups.txt file. The script will produce hic files for viewing in Juicebox as well as scaffold annotation files. Metadata can be attached to the map by passing -i and -g flags with paths to stats and graph files.

ARGUMENTS:
path_to_cprops				Specify cprops file path.
path_to_asm					Specify asm file.
path_to_contig_mnd			Specify path to mnd file describing pairwise Hi-C contacts between assembly contigs.

OPTIONS:
-q mapq						Build map for a specific mapq threshold (default is 1). 
-p true/false				Use GNU Parallel to speed up computation (default is true).
-z zoom						Build map with hardcoded zoom level. By default this is calculated based on the cprops file and applied only to genomes >= 2.1 Gb.

**unprompted**
#-l gap_file				Path to gap bed file - necessary to build contig annotation track.
-m path_to_asm_mnd			Path to mnd already lifted from input to assembly chromosome: used to skip the remapping step.
-n							Skip normalization.
-r							Build for specific resolutions (default is -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000)
-c							Clean up when done (default: no cleanup.)
-i							Ignore mapq suffix.
-h							Shows this help
*****************************************************
"

## Defaults
mapq=1
use_parallel=true
res_string="2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000"

skip_norm=false
clean_up=false
ignore_mapq_suffix=false;
add_metadata=false;

## HANDLE OPTIONS
while getopts "q:p:z:m:l:r:incah" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -q flag was triggered, starting calculations for $OPTARG threshold mapping quality" >&1
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Using the default value MAPQ=${MAPQ[0]}" >&2
        fi
    ;;
	p)	if [ $OPTARG == true ] || [ $OPTARG == false ]; then
    	    echo ":) -p flag was triggered. Running with GNU Parallel support parameter set to $OPTARG." >&1
			use_parallel=$OPTARG
    	else
    		echo ":( Unrecognized value for -p flag. Running with default parameters (-p true)." >&2
    	fi
    ;;
    z)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -z flag was triggered, starting calculations with specified scaling" >&1
            scale=$OPTARG
        else
            echo ":( Wrong syntax for scaling coefficient. Deriving the default scaling to fit the genome into 2.1Gb" >&2
        fi
	;;
	m) if [ -s $OPTARG ]; then 
			echo ":) Skipping remap step and using $OPTARG as premapped input" >&1
			remapped_mnd=$OPTARG
		else
			echo ":( Tentative remapped file not found. Building one as part of the workflow"
		fi		
	;;
#	l)  gap_file=$OPTARG
#	;;
	i)  ignore_mapq_suffix=true;
		echo ":) -i flag was triggered, building mapq without" >&1
	;;
	n)  skip_norm=true
		echo ":) -n flag was triggered, building maps without normalization" >&1
	;;
	c)  clean_up=true
		echo ":) -c flag was triggered, will remove temporary files after completion" >&1
	;;
	r)  re='^[0-9]*(\,?[0-9]*)*$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -r flag was triggered, starting calculations for resolution list: ${OPTARG}" >&1
            res_string=$OPTARG
        else
            echo ":( Wrong syntax for resolution flag. Using the default value pct=5" >&2
        fi
    ;;
    a)  add_metadata=true
		echo ":) -a flag was triggered, will look for juicer metadata files and add if present" >&1
	;;
	*)  echo ":( Illegal options. Exiting."
		echo "$USAGE"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

IFS=',' read -r -a res <<< "$res_string"

## HANDLE ARGUMENTS, TODO: check formats
cprops=$1
asm=$2
mnd=$3
[ $# -eq 3 ] && [ -s ${cprops} ] && [ -s ${asm} ] && [ -s ${mnd} ] || {	echo ":( Not sure how to parse your input or input files not found at intended locations. Exiting!" && echo "$USAGE" && exit 1 ; }

## CHECK DEPENDENCIES
	type java >/dev/null 2>&1 || { echo >&2 ":( Java is not available, please install/add to path Java to run Juicer and Juicebox. Exiting!"; exit 1; }

if [ $use_parallel == true ]; then
	type parallel >/dev/null 2>&1 || { echo >&2 ":( GNU Parallel support is set to true (default) but GNU Parallel is not in the path. Please install GNU Parallel or set -p option to false. Exiting!"; exit 1; }
fi

path_to_scripts=`cd "$( dirname $0)" && pwd`
path_to_lift=$(dirname ${path_to_scripts})"/lift"

juicebox=${path_to_scripts}/"juicebox_tools.sh"

lift_input_mnd_script=${path_to_lift}/lift-input-mnd-to-asm-mnd.awk
lift_input_annotations_script=${path_to_lift}/lift-input-annotations-to-asm-annotations.awk

if [ ! -f $juicebox ] || [ ! -f $lift_input_mnd_script ] || [ ! -f $lift_input_annotations_script ] ; then
    echo ":( Relevant dependency scripts not found. Exiting!" && exit 1
fi


## Handle default zoom if not specified

if [ -z $scale ]; then
	# calculate necessary zoom
        totlength=`awk 'FILENAME==ARGV[1]{len[$2]=$3;len[-$2]=$3;next}{for(i=1;i<=NF;i++){total+=len[$i]}}END{print total}' ${cprops} ${asm}`
	scale=$(( 1 + $totlength / 2100000000 ))
fi

genomeid="`basename $asm .asm`"

if [ -z ${remapped_mnd} ]; then
	## Remap merged_nodups
	echo "...Remapping contact data from the original contig set to assembly"
	if [ $use_parallel == true ]; then
		cmd="parallel --will-cite -a ${mnd} --pipepart -j 80% --block 1G \"awk -v scale=${scale} -f ${lift_input_mnd_script} ${cprops} ${asm} - \" > temp.${genomeid}.asm_mnd.txt"
	else
		cmd="awk -v scale=${scale} -f ${lift_input_mnd_script} ${cprops} ${asm} ${mnd} > temp.${genomeid}.asm_mnd.txt"
	fi

	eval ${cmd}
	remapped_mnd="temp."${genomeid}".asm_mnd.txt"
fi

## Make tracks
echo "...Building track files"

awk 'BEGIN{OFS="\t"; print "chr1", "sx1", "sx2", "chr2", "sy1", "sy2", "color", "Scaffold_ID", "x1", "x2", "y1", "y2"}{print $1, 0, $3, $1, 0, $3, "0,255,0", "+"$1" (+"$2")", 0, $3, 0, $3}' $cprops | awk -v scale=${scale} -f ${lift_input_annotations_script} $cprops $asm - > ${genomeid}_asm.scaffold_track.txt

awk -v scale=${scale} 'BEGIN{OFS="\t"; print "chr1", "sx1", "sx2", "chr2", "sy1", "sy2", "color", "Superscaffold_ID", "x1", "x2", "y1", "y2"; pos+=0}(FILENAME==ARGV[1]){clength[$2]=$3; next}{gsub("-",""); n=split($0,a); c=0; for (i=1; i<=n; i++) {c+=clength[a[i]]}; print "assembly", int(pos/scale), int((pos+c)/scale), "assembly", int(pos/scale), int((pos+c)/scale), "0,0,255", FNR, pos, pos+c, pos, pos+c; pos+=c}' $cprops $asm > ${genomeid}_asm.superscaf_track.txt

# trial single-file for jb4a
awk '{$1=">"$1}1' ${cprops} | cat - ${asm} > ${genomeid}".assembly"

## Build .hic files

echo "...Building the hic file"

[ $mapq -eq 1 ] && ignore_mapq_suffix=true ## lab convention to keep mapq 1 wo suffix, otherwise add suffix in case building multiple maps
[ ${ignore_mapq_suffix} == "true" ] && mapqsuf="" || mapqsuf="_"${mapq}

rLen=${#res[@]}
add_options=$(( res[0]/scale ))
for (( i=1; i<$rLen; i++ ))
do 
    add_options=$add_options","$(( res[$i]/scale ))
done

#[ ${scale} -ne 1 ] && add_options=${res[0]}","${add_options}

add_options="-r "${add_options}

[ "$skip_norm" == "true" ] && add_options=${add_options}" -n"

if [ "$add_metadata" == "true" ]; then
	[ -f inter${mapqsuf}.txt ] && add_options=${add_options}" -s inter${mapqsuf}.txt"
	[ -f inter${mapqsuf}_hists.m ] && add_options=${add_options}" -g inter${mapqsuf}_hists.m"
fi

bash ${juicebox} pre -q ${mapq} ${add_options} ${remapped_mnd} ${genomeid}${mapqsuf}.hic <(echo "assembly	"$((totlength / scale)))


## Cleanup
[ "$clean_up" == "true" ] && rm ${remapped_mnd}

#
#
