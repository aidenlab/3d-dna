#!/bin/bash

#### Description: Wrapper script that applies a list of edits to the original set of scaffolds/contigs to create a new input set, described by new cprops, new mnd and, if requested, new fasta.
#### Usage: apply-edits-prep-for-next-round.sh [ -f path_to_fasta_file ] [ -r revision_label ] [ -p true/false ] <path_to_edit_annotations> <path_to_original_cprops> <path_to_original_mnd>
#### Input: Edits encoded as 2D annotations (0-based, tab-separated, 12 column format), original input cprops and merged_nodups.txt files
#### Dependencies: edit-cprops-according-to-annotations.awk, edit-mnd-according-to-new-cprops.awk, edit-fasta-according-to-new-cprops.awk
#### Options: -r revision label; -f path to original fasta file; -p for GNU Parallel support.
#### Output: New cprops, new mnd file (optional: new fasta file), named as *.revision_label.cprops; *.revision_label.txt
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version date 12/07/2016.

USAGE="
*****************************************************

Usage: ./apply-edits-prep-for-next-round.sh [ -r suffix ] [ -f path_to_fasta_file ] [ -p true/false ] <path_to_edit_annotations> <path_to_original_cprops> <path_to_original_mnd>

ARGUMENTS:
path_to_edit_annotations		Path to 2D annotation files describing regions in original scaffolds/contigs to be labeled as debris.
path_to_original_cprops		Path to cprops describing original input scaffold/contig set
path_to_original_mnd		Path to mnd file describing Hi-C contacts across original intput scaffold/contig set

OPTIONS:
-h				Shows this help
-r	revision_label		Revision label to annotate edited files (default is \"revision\") 
-f	path_to_fasta		Path to fasta file if plan to generate [not sure if needed, not working for now]
-p	true/false		Use GNU Parallel to speed up calculations (default is true)

*****************************************************
"
## SET DEFAULTS
use_parallel="true"
suffix="revision"

## HANDLE OPTIONS
while getopts "r:f:p:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	r)  echo "...-r flag was triggered, output will be labeled as *.${OPTARG}.*" >&1
            suffix="${OPTARG}"
	;;
	f) if [ -s $OPTARG ]; then 
			echo "...-f flag was triggered, will dump edited fasta from $OPTARG" >&1
			orig_fasta=$OPTARG
		else
			echo ":( Original fasta file not found. Continuing without the fasta dump" >&2
		fi
		
	;;
	p)	if [ $OPTARG == true ] || [ $OPTARG == false ]; then
    	    echo "...-p flag was triggered. Running with GNU Parallel support parameter set to $OPTARG." >&1
			use_parallel=$OPTARG
    	else
    		echo ":( Unrecognized value for -p flag. Running with default parameters (-p true)." >&2
    	fi
    ;;
	*)  echo ":( Illegal options. Exiting!"  >&2
		echo "$USAGE"  >&2
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS, TODO: check formats
edits=$1
cprops=$2
mnd=$3

[ $# -eq 3 ] && [ -s ${edits} ] && [ -s ${cprops} ] && [ -s ${orig_mnd} ] || {	echo >&2 ":( Not sure how to parse your input or input files not found at intended locations. Exiting!" && echo "$USAGE" >&2 && exit 1; }

## HANDLE DEPENDENCIES
if [ $use_parallel == true ]; then
	type parallel >/dev/null 2>&1 || { echo >&2 ":( GNU Parallel support is set to true (default) but GNU Parallel is not in the path. Please install GNU Parallel or set -p option to false. Exiting!"; exit 1; }
fi

path_to_scripts=`cd "$( dirname $0)" && pwd`

edit_cprops_script=${path_to_scripts}"/edit-cprops-according-to-annotations.awk"
edit_mnd_script=${path_to_scripts}"/edit-mnd-according-to-new-cprops.awk"
edit_fasta_script=${path_to_scripts}"/edit-fasta-according-to-new-cprops.awk"

if [ ! -f ${edit_cprops_script} ] || [ ! -f ${edit_mnd_script} ] || [ ! -f ${edit_fasta_script} ]; then
    echo >&2 ":( Relevant dependency scripts not found. Exiting!" && exit 1
fi

## MAIN

## edit cprops
echo "...applying edits to cprops file" >&1
filename=$(basename "$cprops")
extension="${filename##*.}"
filename="${filename%.*}"
awk -v label1=":::fragment_" -v label2=":::debris" -f ${edit_cprops_script} ${edits} ${cprops} > ${filename}.${suffix}.${extension}
new_cprops="${filename}.${suffix}.${extension}"

## edit mnd
echo "...applying edits to mnd file" >&1
filename=$(basename "$mnd")
extension="${filename##*.}"
filename="${filename%.*}"
if [ ${use_parallel} == "true" ]; then
	parallel -a ${mnd} --pipepart --will-cite --jobs 80% --block 1G "awk -v label1=\":::fragment_\" -v label2=\":::debris\" -f ${edit_mnd_script} ${new_cprops} - " > ${filename}.${suffix}.${extension}
else	
	awk -v label1=":::fragment_" -v label2=":::debris" -f ${edit_mnd_script} ${new_cprops} ${mnd} > ${filename}.${suffix}.${extension}
fi

## edit fasta
if [ ! -z ${fasta} ] && [ -s ${fasta} ]; then
	echo "...applying edits to input sequence file" >&1
	filename=$(basename "$fasta")
	extension="${filename##*.}"
	filename="${filename%.*}"
	awk -v label1=":::fragment_" -v label2=":::debris" -f ${edit_fasta_script} ${new_cprops} ${fasta} > ${filename}.${suffix}.${extension}
fi

