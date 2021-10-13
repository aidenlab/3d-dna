#!/bin/bash
# compact version of the script to extract Hi-C reads that overlap SNPs.
# DEPENDENCY: Gnu parallel wrapper for extract-edges awk script (optional but recommended)
# Written by: OD (olga.dudchenko@bcm.edu)
# Original version: 2014, latest mods 210122.

USAGE="
**********************

Script to extract from the merged_nodups file (long format) Hi-C contacts that overlap SNPs.

	./extract-SNP-reads-from-mnd-file.sh [-q mapq] [-n] [-d] [-a] [-c] [-h] [-o output_file_name] <path_to_psf_file> <path_to_mnd_txt>

Input: psf file (\">CHR POS VAR1 VAR2 ID ...\") and mnd file.

Options:
	-q  <mapq>:
                data from reads with mapping quality less than mapq will be ignored. Default: 1.
    -n
                keeps 'native' coordinates in the mnd file. Default: coordinates are lifted for the benefit of the SNP-based contact maps.
    -d
                keeps 'dangling' edges (read pairs in which only one read overlaps with a SNP or SNPs.
    -a
                amplifies edges that correspond to SNPs that 'live' on the same read. This is to increase the weight of the relevant edges (typically should already be phased by the SNP caller).
    -c
                reports edges combinatorially in cases when single reads overlap multiple SNPs. 
	-o  <outfile_file_name>
                filename to dump data. Default: add snp suffix to name of mnd file, e.g. mnd.txt -> mnd.snp.txt
    -h
                shows this help

************************
"
## DEFAULTS
mapq=1
keep_native_coordinates=0
keep_dangling_edges=0
amplify_neighboring_snps=0
report_combinatorial_edges=0

## HANDLE OPTIONS
while getopts "q:o:ndach" opt; do
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
    n)
        keep_native_coordinates=1;
    ;;
    d)
        keep_dangling_edges=1;
    ;;
    a)
        amplify_neighboring_snps=1
    ;;
    c)
        report_combinatorial_edges=1
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
single_thread_script=${pipeline}/extract-SNP-reads-from-mnd-file.awk

[ $parallel == "true" ] && parallel -a ${mnd_file} --pipepart --will-cite --jobs 50% --block 1G "awk -v mapq=${mapq} -v keep_native_coordinates=${keep_native_coordinates} -v keep_dangling_edges=${keep_dangling_edges} -v amplify_neighboring_snps=${amplify_neighboring_snps} -v report_combinatorial_edges=${report_combinatorial_edges} -f ${single_thread_script} ${psf_file} - " > ${outfile} || awk -v mapq=${mapq} -v keep_native_coordinates=${keep_native_coordinates} -v keep_dangling_edges=${keep_dangling_edges} -v amplify_neighboring_snps=${amplify_neighboring_snps} -v report_combinatorial_edges=${report_combinatorial_edges} -f ${single_thread_script} ${psf_file} ${mnd_file} > ${outfile}
