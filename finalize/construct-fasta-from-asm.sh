#!/bin/bash
## Wrapper script to generate a fasta from internal datatypes: cprops and asm
## Input: internally consistent cprops, asm and fasta
## TODO: parallelize for speed
## Prints into STDOUT

USAGE="
*****************************************************
USAGE: construct-fasta-from-asm.sh <path-to-cprops> <path-to-asm> <path-to-fasta>
*****************************************************
"

## HANDLE OPTIONS

label="HiC_scaffold"

while getopts "l:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
    l) label=$OPTARG
    ;;
	*) echo "$USAGE" >&2
		exit 1
	;;
esac
done
shift $(( OPTIND-1 ))


# DEPENDENCIES:
pipeline=`cd "$( dirname $0)" && cd .. && pwd`


## HANDLE ARGUMENTS
[ -z $1 ] || [ -z $2 ] || [ -z $3 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

cprops=$1
asm=$2
fasta=$3

awk -f ${pipeline}/utils/index-fasta.awk $fasta > tmp.index.txt
declare -A index
while read name startbyte skip
do
	index[$name]=$startbyte
done < <(awk 'FILENAME==ARGV[1]{cname[$1]=$2;next}{$1=cname[$1]; print}' $cprops tmp.index.txt)


# read asm line by line
scaffold_counter=1
while read -r line
do
	echo ">"${label}"_"${scaffold_counter}
		
	for contig in $line
	do
		## might use more annotations
		if [[ $contig == -* ]]; then
			contig=$(echo $contig | cut -d'-' -f 2)
			tail -c +${index[${contig}]} ${fasta} | awk '$0~/>/{exit}1' | awk -f ${pipeline}/utils/reverse-fasta.awk -
		else
			tail -c +${index[${contig}]} ${fasta} | awk '$0~/>/{exit}1'
		fi
	done
	let scaffold_counter=scaffold_counter+1  		
done < $asm

rm tmp.index.txt
