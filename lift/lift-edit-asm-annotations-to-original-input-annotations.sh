#!/bin/bash

#### Description: Does a series of liftovers to map annotations from an edited assembly to original scaffold/contig input.
#### Usage: lift-edit-asm-annotations-to-original-input-annotations.sh <path_to_original_cprops> <path_to_edited_cprops> <path_to_edited_asm> <path_to_edited_annotation_file>
#### Input: Current and original cprops, current asm file.
#### Optional input: list of previous edits to original contigs/scaffolds.
#### Output: Stdout of original contig/scaffold edits as a 2D annotation file (0-based, 12-column).
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version date 01/11/2017.

USAGE="
*****************************************************
Translates list of assembly mismatch positions into mismatches in original contigs/scaffolds (and merges them with previous annotations, optional).

./lift-edit-asm-annotations-to-original-input-annotations.sh [-s scale] <path_to_original_cprops> <path_to_edited_cprops> <path_to_edited_asm> <path_to_edited_annotation_file>

ARGUMENTS:
path_to_original_cprops		Path to original input cprops file.
path_to_edited_cprops		Path to cprops file of the assembly made from edited input
path_to_edited_asm		Path to asm file of the assembly made from edited input
path_to_edited_annotation_file		Path to a 2D annotation file to be lifted from current edited assembly to original input scaffolds/contigs


OPTIONS:
-h				Shows this help
#-s	scale		Scaling coefficient to stretch or squeeze the annotations [currently hard-set to default 1].
*****************************************************
"

## SET DEFAULTS
scale=1

## HANDLE OPTIONS
while getopts "s:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	s)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo "...-s flag was triggered. Will scale the output according to $OPTARG scaling coefficient" >&2
            scale="${suffix}${OPTARG}"
        else
            echo ":( Wrong syntax for scaling coefficient. Continuing without scaling" >&2
        fi
	;;
	*)  echo ":( Illegal options. Exiting."
		echo "$USAGE"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS, TODO: check formats
orig_cprops=$1
edit_cprops=$2
edit_asm=$3
annotations=$4

[ $# -eq 4 ] && [ -s ${orig_cprops} ] && [ -s ${edit_cprops} ] && [ -s ${edit_asm} ] && [ -s ${annotations} ] || {	echo >&2 ":( Not sure how to parse your input or input files not found at intended locations. Exiting!" && echo "$USAGE" >&2 && exit 1; }

## HANDLE DEPENDENCIES
path_to_scripts=`cd "$( dirname $0)" && pwd`
path_to_lift=$(dirname ${path_to_scripts})"/lift"

lift_asm_annotations_script=${path_to_lift}"/lift-asm-annotations-to-input-annotations.awk"
lift_input_annotations_script=${path_to_lift}"/lift-input-annotations-to-asm-annotations.awk"

## MAIN
head -n 1 ${annotations} && awk '{print $2}' ${orig_cprops} | awk -v scale=${scale} -f ${lift_asm_annotations_script} ${orig_cprops} - <(awk '{print $2}' ${edit_cprops} | awk -v scale=${scale} -f ${lift_input_annotations_script} ${edit_cprops} - <(awk -f ${lift_asm_annotations_script} ${edit_cprops} ${edit_asm} ${annotations})) | tail -n +2 | sort -k 1,1 -k 2,2n
