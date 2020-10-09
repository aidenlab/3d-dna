#!/bin/bash
# Helper script to rapidly generate fasta subset using line index
# INPUT: comma-separated list of sequence names to be printed, path to master fasta file
# OUTPUT: stdout in standard fasta format
# OPTIONAL: index file (-i <index>)
# Written by OD, 190223

USAGE="
*****************************
	./print-fasta-subset.sh [ -i <path_to_index_file_or_named_pipe> ] <comma_separated_list_of_sequence_names_to_print> <path_to_master_fasta>
*****************************
"
while getopts "i:h" opt; do
case $opt in
    i) if [ -f ${OPTARG} ] || [ -p ${OPTARG} ]; then
            index=$OPTARG
        fi
    ;;
    h) echo "$USAGE" >&1
        exit 0
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

list=$1
fasta=$2

## HANDLE ARGUMENTS
[ -z $1 ] || [ -z $2 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

IFS=',' read -r -a array <<< "$list"
for element in "${array[@]}"
do
	if [ -z ${index} ]; then
		# could do faster but perhaps worthwhile keeping order
		awk -v sequence="$element" '$0~/^>/{test=0}$1==">"sequence{test=1}test' ${fasta}
	else
		IFS=' ' read skip firstline lastline <<< $(awk -v sequence="$element" '$1==sequence' ${index})
		cmd="sed -n '${firstline},${lastline}p;$((lastline+1))q' ${fasta}"
		eval $cmd
	fi
done

