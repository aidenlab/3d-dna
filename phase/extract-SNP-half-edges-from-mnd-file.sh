#!/bin/bash
# parallelization for script to extract Hi-C contacts that overlap SNPs.
# DEPENDENCY: Gnu parallel wrapper for extract-edges awk script (optional but recommended)
# Written by: OD (olga.dudchenko@bcm.edu)
# Original version: 2014, latest mods 2019.

USAGE="
**********************

This is a parallel verion of the script to extract Hi-C contacts that overlap SNPs.

	./extract-SNP-edges-from-mnd-file.sh [-q mapq] [-o output_file_name] <path_to_sorted_psf_file> <path_to_mnd_txt>

Input: psf file (\">CHR POS VAR1 VAR2 ...\") and mnd file.

Options:
	-q <mapq>: data from reads with mapping quality less than mapq will be ignored. Default: 1.
	-o <outfile_file_name>: filename to dump data. Default: add snp suffix to name of mnd file, e.g. mnd.txt -> mnd.snp.txt

**********************
"
## DEFAULTS
mapq=1

## HANDLE OPTIONS
while getopts "q:o:h" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
    q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo "-q flag was triggered, setting mapping quality threshold to $OPTARG." >&1
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality parameter. Exiting!" >&2
            exit 1
        fi
    ;;
    o)
    	outfile=$OPTARG
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS

psf_file=$1
mnd_file=$2

([ -z ${psf_file} ] || [ -z ${mnd_file} ]) && (":( Some input is missing. Exiting!" >&2 && exit 1)

[ -z ${outfile} ] && outfile=`basename ${mnd_file} .txt`".snp.txt"

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo >&2 ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!"

pipeline=`cd "$( dirname $0)" && pwd`
single_thread_script=${pipeline}/extract-SNP-half-edges-from-mnd-file.awk

[ $parallel == "true" ] && parallel -a ${mnd_file} --pipepart --will-cite --jobs 80% --block 1G "awk -v mapq=${mapq} -f ${single_thread_script} ${psf_file} - " > ${outfile} || awk -v mapq=${mapq} -f ${single_thread_script} ${psf_file} ${mnd_file} > ${outfile}
