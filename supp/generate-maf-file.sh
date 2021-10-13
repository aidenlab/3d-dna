#!/bin/bash
## This is a helper script to generate whole-genome alignment in maf format 
## Usage: ./generate-maf-file.sh [-s split_size] <target_fasta> <query_fasta>
## Written by: Olga Dudchenko

USAGE=" ./generate-maf-file.sh [ -a alignment_score_threshold ] <target_fasta> <query_fasta> 
	number_of_threads	number of threads for parallelization

OPTIONS:
	prompted
	-a alignment_threshold	passed on to lastz as hspthres parameter [default: 50000]

	unprompted:
	-s split_size	parallelization parameter for lenght to be handled by single threads [default: 300000000]
"

# HANDLE DEPENDENCIES

## 3D-DNA
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## LASTZ dependency
lastz="false"
if hash lastz 2>/dev/null; then
	lastz="true"
fi
[[ $lastz == "false" ]] && echo >&2 ":( LASTZ not installed or not in the path. The merge section of the pipeline will not run. Exiting!" && exit 1

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Not supporting this script wo parallel. Exiting!" >&2 && exit 1

# DEFAULTS
splitsize=300000000
hspthresh=50000

# HANDLE OPTIONS

re='^[0-9]+$'

while getopts "s:a:h" opt; do
case $opt in
	a) 	if [[ $OPTARG =~ $re ]]; then
			echo " -a flag was triggered, alignment threshold set to $OPTARG." >&1
			hspthresh=$OPTARG
		else
			echo ":( Wrong syntax for flag. Exiting!" >&2
			exit 1
		fi
	;;
	s) 	if [[ $OPTARG =~ $re ]]; then
			echo " -s flag was triggered, fasta length per thread set to $OPTARG." >&1
			splitsize=$OPTARG
		else
			echo ":( Wrong syntax for flag. Exiting!" >&2
			exit 1
		fi
	;;
	h) echo "$USAGE"
	exit 0
	;;
	*) echo "$USAGE" >&2
	exit 1
	;;
esac
done
shift $(( OPTIND-1 ))

# HANDLE ARGUMENTS

target_fasta=$1
query_fasta=$2

([ -f ${target_fasta} ] && [ -f ${query_fasta} ]) || (echo ":( Expected arguments not found. Exiting!" && exit 1)

targetid=$(basename "$target_fasta" .fasta)
targetid=$(basename "$targetid" .fna)
targetid=$(basename "$targetid" .fa)

queryid=$(basename "$query_fasta" .fasta)
queryid=$(basename "$queryid" .fna)
queryid=$(basename "$queryid" .fa)

[ -z ${target_chrom} ] && ln -sf ${target_fasta} ${targetid}.tmp.fasta || awk -f ${asm_pipeline}/finalize/make-fasta-subset.awk ${target_chrom} ${target_fasta} > ${targetid}.tmp.fasta
[ -z ${query_chrom} ] && ln -sf ${query_fasta} ${queryid}.tmp.fasta || awk -f ${asm_pipeline}/finalize/make-fasta-subset.awk ${query_chrom} ${query_fasta} > ${queryid}.tmp.fasta

# MAIN

# 1) to speed things up split up the target genome into chunks

n=`awk -v ignore_sorting=1 -f ${pipeline}/utils/generate-assembly-file-from-fasta.awk ${targetid}.tmp.fasta | awk '$0~/^>/{print substr($0,2)}' | awk -v splitsize=${splitsize} 'BEGIN{counter++}{if(c>splitsize){counter++; c=0}; print $1 >"split."counter".txt"; c+=$NF;}END{print counter}' -`

echo $n

seq 1 $n | parallel "awk -f ${pipeline}/finalize/make-fasta-subset.awk split.{}.txt ${targetid}.tmp.fasta > ${targetid}.{}.fasta"
rm ${targetid}.tmp.fasta

# 2) actual alignment

seq 1 $n | parallel "lastz ${targetid}.{}.fasta[multiple] ${queryid}.tmp.fasta --notransition --step=20 --nogapped --format=maf --ambiguous=iupac --hspthresh=${hspthresh}" > ${queryid}_to_${targetid}.maf

# 3) cleanup
seq 1 $n | parallel "rm split.{}.txt"
seq 1 $n | parallel "rm ${targetid}.{}.fasta"
rm ${queryid}.tmp.fasta
