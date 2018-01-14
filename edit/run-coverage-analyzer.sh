#!/bin/bash
#### Description: Wrapper script to annotate likely repeats. Analyzes normalization vector at a specified resolution (default 25kb).
#### Usage: run-mismatch-detector.sh -w <bin_size> <path-to-hic-file>
#### Dependencies: Juicebox_tools
#### Input: Juicebox hic file.
#### Output: "Wide" bed file highlighting likely repeat regions [repeat_wide.bed; maybe later: repeat_narrow.bed]. Additional output generated as part of this wrapper script includes repeat_score_wide.wig (repeat_score_narrow.wig) track files.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version date 12/19/2016.

USAGE="
*****************************************************
This is a wrapper for a fragment of Hi-C misassembly detection pipeline, version date: Dec 3, 2016. This fragment concerns with generating a mismatch annotation file that will later be overlaid with scaffold boundaries to excise regions spanning misassemblies.

Usage: ./run-mismatch-detector.sh [-h] [-p percentile] [-b bin_size_aka_resolution] [-d depletion_region_size] path_to_hic_file

ARGUMENTS:
path_to_hic_file     	Path to Juicebox .hic file of the current assembly.

OPTIONS:
-h			Shows this help
-w wide_res	Sets resolution for the first-pass search of repeats (default is 25000 bp)
-t thr_cov	Threshold coverage [default is 2]

Unprompted:
...

Uses juicebox_tools.sh that should be in the same folder as the wrapper script.
*****************************************************
"

## Set defaults
bin_size=25000			# default bin size to do a first-pass search for repeats
thr_cov=2

## HANDLE OPTIONS

while getopts "hw:t:" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
    w)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -w flag was triggered, performing cursory search for repeat at $OPTARG resolution" >&1
            bin_size=$OPTARG
        else
            echo ":( Wrong syntax for bin size. Using the default value 25000" >&2
        fi
    ;;
    t)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -t flag was triggered, flagging regions with coverage higher than $OPTARG" >&1
            thr_cov=$OPTARG
        else
            echo ":( Wrong syntax for threshold coverage. Using the default value 2" >&2
        fi
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS: TODO check file format
if [ $# -lt 1 ]; then
    echo ":( Required arguments not found. Please double-check your input!!" >&2
    echo "$USAGE" >&2
    exit 1
fi

hic_file=$1

## CHECK DEPENDENCIES
	type java >/dev/null 2>&1 || { echo >&2 ":( Java is not available, please install/add to path Java to run Juicer and Juicebox. Exiting!"; exit 1; }

path_to_scripts=`cd "$( dirname $0)" && pwd`
path_to_vis=$(dirname ${path_to_scripts})"/visualize"

juicebox=${path_to_vis}/"juicebox_tools.sh"

## DUMP NV FOR ANALYSIS, ANNOTATE ALL >=2. TODO: check that the matrix is not too sparse for chosen resolution
echo "...Dumping ${bin_size} resolution coverage track"

bash ${juicebox} dump norm VC ${hic_file} assembly BP ${bin_size} "coverage_wide.wig"
[ $? -ne 0 ] && echo >&2 ":( Juicebox coverage dump is empty! Perhaps something is wrong with the hic file or the requested resolution is too high. Exiting!" && exit 1

echo "fixedStep chrom=assembly start=1 step="${bin_size}" span="${bin_size} | cat - "coverage_wide.wig" > temp && mv temp "coverage_wide.wig"

## NOTE: if(start) was added after 3d-dna release. Caused bug when nothing was annotated due to coverage.
awk -v thr_cov=${thr_cov} -v bin_size=${bin_size} 'NR>1&&$0>=thr_cov{print (NR-2)*bin_size, (NR-1)*bin_size}' coverage_wide.wig | awk 'BEGIN{OFS="\t"}NR==1{start=$1; end=$2;next}$1==end{end=$2;next}{print "assembly", start, end; start=$1; end=$2}END{if(start){print "assembly", start, end}}' > repeats_wide.bed

## TODO: maybe downstread add filtering by looking for nearby mismatches to get rid of false positives
